# Random linear physical networks

Code accompanying the paper<br>
**Impact of physicality on network structure**<br>
Márton Pósfai, Balázs Szegedy, Iva Bačić, Luka Blagojević, Miklós Abért,
János Kertész, László Lovász, Albert-László Barabási
[arXiv:2211.13265](https://arxiv.org/abs/2211.13265)

The repository contains three elements:

1. `randLPN.py`: a minimal package implemented in pure python, which allows to generate random linear physical networks.
2. `randLPN_CPP.cpp`: a python extension written in c++, also allowing to generate random linear physical networks. Faster than the python implementation. To build the package run: `python setup.py build`
4. `figures/`: a collection of Jupyter notebooks that reproduce the four figures in the main text of arXiv:2211.13265.
