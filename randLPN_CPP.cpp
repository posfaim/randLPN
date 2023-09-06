//Python module to generate linear physical networks

#define PY_SSIZE_T_CLEAN
#include <Python.h>

#include <iostream>
#include <chrono>
#include <random>

#include <CGAL/Simple_cartesian.h>

#define SMALL_NUMBER 1e-14
#define MAX_NODE_TRIAL 100000
#define M_PI           3.14159265358979323846  /* pi */

typedef CGAL::Simple_cartesian<double> Kernel;
typedef Kernel::Point_3 PT;
typedef Kernel::Segment_3 SEG;


// construct a trivial random generator engine from a time-based seed:
unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
std::default_random_engine generator(seed);

// 
template <class GEOM1, class GEOM2>
bool DistTest(double squared_dist, GEOM1& geom, std::vector<GEOM2> &geoms) {
    for (GEOM2 geom2: geoms){
        double sd = CGAL::squared_distance(geom, geom2);
        if ((sd>SMALL_NUMBER) && (sd<squared_dist)){ // questionable practice
            return false;
        }
    }
    return true;
}

// Physical ER
int RandPointCloud(int N, double dist, std::vector<PT> &points){
    int max_node_trial = MAX_NODE_TRIAL;
    double squared_dist=dist*dist;
    PT point;
    
    std::uniform_real_distribution<double> distribution(0.0,1.0);
    
    int n;
    for(n=0;n<N;n++){
        int trial = 0;
        do{
            point = PT(distribution(generator),distribution(generator),distribution(generator));
            trial++;
        }while( (!DistTest<PT,PT>(squared_dist, point, points)) && (trial<max_node_trial));

        if (trial < max_node_trial){
            points.push_back(point);
        }else{
            n--;
            break;
        }
    }
    return n;
}


int LinPhysER_C(int N, int L, double dist, std::vector<PT> &points, std::vector<std::array<int,2>> &vws, std::vector<int> &t_add){
    int n = RandPointCloud(N, dist, points);
    double squared_dist=dist*dist;
    
    //edges are attempted to be added in random order
    std::vector<std::array<int,2>> order;
    order.reserve(n*(n-1)/2);
    for(int v=0;v<n;v++)for(int w=v+1;w<n;w++) order.push_back({v,w});
    std::shuffle(std::begin(order), std::end(order), generator);
    
    //let's get ready to add those edges
    int l = 0, t = 0;
    std::vector<SEG> segs;
    segs.reserve(L);
    vws.reserve(L);
    t_add.reserve(L);
    
    //let's add those edges
    for (auto vw: order){
        SEG seg(points[vw[0]],points[vw[1]]);
        if ( DistTest<SEG,PT>(squared_dist, seg, points) && DistTest<SEG,SEG>(squared_dist, seg, segs)){
            segs.push_back(seg);
            vws.push_back({vw[0],vw[1]});
            t_add.push_back(t);
            l++;
            if (l==L) break;
        }
        t++;
    }
    
    return n;    
}

//only considering link-link interaction
//nodes are still placed at least lambda (=dist) away from each other
int LinPhysER_linkonly_C(int N, int L, double dist, std::vector<PT> &points, std::vector<std::array<int,2>> &vws, std::vector<int> &t_add){
    double squared_dist=dist*dist;
    
    std::uniform_real_distribution<double> distribution(0.0,1.0);
    for(int n=0;n<N;n++)
        points.push_back(PT(distribution(generator),distribution(generator),distribution(generator)));
    
    //edges are attempted to be added in random order
    std::vector<std::array<int,2>> order;
    order.reserve(N*(N-1)/2);
    for(int v=0;v<N;v++)for(int w=v+1;w<N;w++) order.push_back({v,w});
    std::shuffle(std::begin(order), std::end(order), generator);
    
    //let's get ready to add those edges
    int l = 0, t = 0;
    std::vector<SEG> segs;
    segs.reserve(L);
    vws.reserve(L);
    t_add.reserve(L);
    
    //let's add those edges
    for (auto vw: order){
        SEG seg(points[vw[0]],points[vw[1]]);
        if (DistTest<SEG,SEG>(squared_dist, seg, segs)){
            segs.push_back(seg);
            vws.push_back({vw[0],vw[1]});
            t_add.push_back(t);
            l++;
            if (l==L) break;
        }
        t++;
    }
    
    return N;    
}

// distance for all pairs of segments upto a threshold
void AllToAllSegmentSquaredDist_C(int N, double maxsd, std::vector<PT> &points, std::vector<std::array<int,2>> &vws,std::vector<std::array<int,2>> &metaedges, std::vector<double> &squareddists){
    //pick random points
    std::uniform_real_distribution<double> distribution(0.0,1.0);
    for(int n=0;n<N;n++)
        points.push_back(PT(distribution(generator),distribution(generator),distribution(generator)));
    
    // create segment objects
    std::vector<SEG> segs;
    segs.reserve(N*(N-1)/2);
    vws.reserve(N*(N-1)/2);
    
    for(int v=0;v<N;v++) for(int w=v+1;w<N;w++){
        segs.push_back(SEG(points[v],points[w]));
        vws.push_back({v,w});
    }
    
    double sd;
    for(int s=0;s<int(vws.size());s++)for(int t=s+1;t<int(vws.size());t++){
        if(vws[s][0]!=vws[t][0] && vws[s][0]!=vws[t][1] && vws[s][1]!=vws[t][0] && vws[s][1]!=vws[t][1]){
            sd = CGAL::squared_distance(segs[s], segs[t]);
            if (sd<maxsd){
                squareddists.push_back(sd);
                metaedges.push_back({s,t});
            }
        }
    }
    
    return;
}


//****************************
//Exposing functions to python

static PyObject* LinPhysER(PyObject* self, PyObject* args){
    int N, L;
    double d;
    
    
    if (!PyArg_ParseTuple(args, "iid", &N, &L, &d))
        return NULL;       
    
    std::vector<PT> points;
    points.reserve(N);
    std::vector<std::array<int,2>> vws;
    std::vector<int> t_add;
        
    LinPhysER_C(N, L, d, points, vws, t_add);

    PyObject* python_points = PyList_New(points.size());
    for (int v = 0; v<int(points.size()); v++)
        PyList_SetItem(python_points, v, Py_BuildValue("(fff)",points[v].x(),points[v].y(),points[v].z()));
    
    PyObject* python_vws = PyList_New(vws.size());
    for (int e = 0; e<int(vws.size()); e++)
        PyList_SetItem(python_vws, e, Py_BuildValue("(ii)", vws[e][0], vws[e][1]));
    
    PyObject* python_t_add = PyList_New(t_add.size());
    for (int e = 0; e<int(t_add.size()); e++)
        PyList_SetItem(python_t_add, e,  PyLong_FromLong(t_add[e]));
    
    
    //return Py_BuildValue("(OO)",python_points,python_vws);
    PyObject* return_tuple = Py_BuildValue("(OOO)",python_points,python_vws, python_t_add);
    
    Py_DECREF(python_points);
    Py_DECREF(python_vws);
    Py_DECREF(python_t_add);
    
    return return_tuple;
}

static PyObject* LinPhysER_linkonly(PyObject* self, PyObject* args){
    int N, L;
    double d;
    
    if (!PyArg_ParseTuple(args, "iid", &N, &L, &d))
        return NULL;       
    
    std::vector<PT> points;
    points.reserve(N);
    std::vector<std::array<int,2>> vws;
    std::vector<int> t_add;
    
    LinPhysER_linkonly_C(N, L, d, points, vws, t_add);

    PyObject* python_points = PyList_New(points.size());
    for (int v = 0; v<int(points.size()); v++)
        PyList_SetItem(python_points, v, Py_BuildValue("(fff)",points[v].x(),points[v].y(),points[v].z()));
    
    PyObject* python_vws = PyList_New(vws.size());
    for (int e = 0; e<int(vws.size()); e++)
        PyList_SetItem(python_vws, e, Py_BuildValue("(ii)", vws[e][0], vws[e][1]));
    
    
    PyObject* python_t_add = PyList_New(t_add.size());
    for (int e = 0; e<int(t_add.size()); e++)
        PyList_SetItem(python_t_add, e,  PyLong_FromLong(t_add[e]));
    
    
    //return Py_BuildValue("(OO)",python_points,python_vws);
    PyObject* return_tuple = Py_BuildValue("(OOO)",python_points,python_vws, python_t_add);
    
    Py_DECREF(python_points);
    Py_DECREF(python_vws);
    Py_DECREF(python_t_add);
    
    return return_tuple;
}

static PyObject* AllToAllSegmentSquaredDist(PyObject* self, PyObject* args){
    int N;
    double maxsd;
    
    std::cout<<"hm1"<<std::endl;
    
    if (!PyArg_ParseTuple(args, "id", &N, &maxsd))
        return NULL;
    
    std::vector<PT> points;
    points.reserve(N);
    std::vector<std::array<int,2>> vws;
    std::vector<std::array<int,2>> metaedges;
    std::vector<double> squareddists;
    
    AllToAllSegmentSquaredDist_C(N, maxsd, points, vws, metaedges, squareddists);
    
    PyObject* python_points = PyList_New(points.size());
    for (int v = 0; v<int(points.size()); v++)
        PyList_SetItem(python_points, v, Py_BuildValue("(fff)",points[v].x(),points[v].y(),points[v].z()));
    
    PyObject* python_vws = PyList_New(vws.size());
    for (int e = 0; e<int(vws.size()); e++)
        PyList_SetItem(python_vws, e, Py_BuildValue("(ii)", vws[e][0], vws[e][1]));
    
    PyObject* python_metaedges = PyList_New(metaedges.size());
    for (int e = 0; e<int(metaedges.size()); e++)
        PyList_SetItem(python_metaedges, e, Py_BuildValue("(ii)", metaedges[e][0], metaedges[e][1]));
    
    PyObject* python_squareddists = PyList_New(squareddists.size());
    for (int x = 0; x<int(squareddists.size()); x++)
        PyList_SetItem(python_squareddists, x, Py_BuildValue("f",squareddists[x]));
    
    PyObject* return_tuple = Py_BuildValue("(OOOO)",python_points,python_vws,python_metaedges,python_squareddists);
    
    Py_DECREF(python_points);
    Py_DECREF(python_vws);
    Py_DECREF(python_metaedges);
    Py_DECREF(python_squareddists);
    
    return return_tuple;
}


static PyMethodDef mainMethods[] = {
    {"LinPhysER",LinPhysER, METH_VARARGS,"Generate linear physical ER"},
    {"LinPhysER_linkonly",LinPhysER_linkonly, METH_VARARGS,"Generate linear physical ER, link-link interaction only"},
    {"AllToAllSegmentSquaredDist",AllToAllSegmentSquaredDist, METH_VARARGS,"Drop N random points and calculate the distance between all segment pairs"},
    {NULL,NULL,0,NULL}
};

static PyModuleDef randLPN_CPP = {
    PyModuleDef_HEAD_INIT,
    "randLPN_CPP","randLPN_CPP",
    -1,
    mainMethods
};

PyMODINIT_FUNC PyInit_randLPN_CPP(void){
    return PyModule_Create(&randLPN_CPP);
}

