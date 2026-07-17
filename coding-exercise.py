"""DIS in real life.

This is a modified version of the exercises developed for the Advanced Artificial Intelligence for precision High Energy Physics school 2023 in Como - see https://github.com/NNPDF/como-2023.
This file derived from a Jupyter notebook so proceed in reading in steps of `In[$]`.
You may also want to comment the later cells out until you reach them.

The script requires a recent version of Python3 and the following additional packages:
- numpy
- scipy
- pandas
- matplotlib
- lhapdf
and in addition the PDF "NNPDF40_nnlo_as_01180".

The LHAPDF library and the actual PDF are available from https://www.lhapdf.org/.
"""

# In[1]:
from collections.abc import Callable
from dataclasses import dataclass
import warnings
import lhapdf
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.integrate import quad, IntegrationWarning
from scipy.special import zeta

# We find at NLO accuracy
# $$
#   \frac 1 x F_2(x,Q^2) = \sum_a (C_a \otimes f_a^p)(x,Q^2) = \sum_a \int\limits_x^1 \frac{dz}{z} C_a(z,Q^2) f_a^p(x/z,Q^2)  = \sum_a\sum_{n=0} \left(\frac{\alpha_s(Q^2)}{4\pi}\right)^n \int\limits_x^1 \frac{dz}{z} C_a^{(n)}(z) f_a^p(x/z,Q^2)
# $$
# with $C_a(z) = e_a^2 c_q(z)$ for the quark channels, i.e. $a=u,\bar u,d,\bar d,s,\bar s,\ldots$,
# and $C_g(z) = e_{tot}^2 c_g(z)$ for the gluon channel, with the total squared charge $e_{tot}^2$.
# The raw coefficient functions are given by
# $$ c_q^{(0)}(z) = \delta(1-z) \quad c_g^{(0)}(z) = 0 $$
# $$ c_q^{(1)}(1) = C_F \left( 4 \left(\frac{\ln(1-z)}{1-z}\right)_+ - 3 \left(\frac{1}{1-z}\right)_+ - (6+4\zeta(2))\delta(1-z) -2(1+z)\ln((1-z)/z) -4 \frac{\ln(z)}{1-z} +6 +4z \right)  $$
# $$ c_g^{(1)}(z) = \left((2-4z(1-z)\ln((1-z)/z)) -2 +16z(1-z)\right) $$
#
# Note that with respect to the lecture notes we have reshuffled the factor $4\pi$ and that the quoted coefficient function actually refer to $1/x F_2(x,Q^2)$.

# When looking at our master equation from above, we immediately face the problem that distributions, such as Dirac delta distributions and plus distributions, are abstract mathematical objects. Recall the definition of plus distributions
# $$\int\limits_0^1 dz g(z) (c(z))_+ = \int\limits_0^1 dz (g(z) - g(1)) c(z) \quad$$
# and you note that this form does not coincide with the factorization formula (note the integration limits!) - we thus need to find a different representation.
#
# Distributions are a common feature of NLO calculation as they have a one-to-one mapping with actual physics features such as soft and collinear divergencies. How to deal with these singularities is discussed under the heading of subtraction schemes.

# Here, a suitable representation is the "RSL scheme" where we divide any coefficient function $c$ in a **R**egular, **S**ingular, and **L**ocal part defined by their behavior under the convolution integral
# $$ (c\otimes f)(x) = \int\limits_x^1 \frac{dz}{z} f(x/z) c^R(z) +  \int\limits_x^1 dz\left(\frac{f(x/z)}{z} - f(x)\right)  c^S(z) + f(x)c^L(x)$$
# where now $c(z)$ may or may not contain distributions.


# ### NLO DIS coefficient functions

# In[2]:


# let's fix some global constants first
CF = 4.0 / 3.0  # Casimir constant in fundamental color representation
e2u = 4.0 / 9.0  # squared electric charge of up-like quarks
e2d = 1.0 / 9.0  # squared electric charge of down-like quarks
nf = 5  # number of active (light) flavors


# Most of the numbers should look familiar to you - we will comment in the next section about the number of light flavors.
#
# Now, we can give the equivalent RSL representation of the NLO DIS coefficient functions:

# In[3]:


# Let's define a helper class to collect all ingredients to the RSL representation: three functions
@dataclass(frozen=True)
class RSL:
    regular: Callable[[float], float]
    singular: Callable[[float], float]
    local: Callable[[float], float]


# we can now spell out all coeffiecient functions
# LO quark
def cq_LO_R(z: float) -> float:
    return 0.0


def cq_LO_S(z: float) -> float:
    return 0.0


def cq_LO_L(x: float) -> float:
    return 1.0


cq_LO = RSL(regular=cq_LO_R, singular=cq_LO_S, local=cq_LO_L)


# NLO quark
def cq_NLO_R(z: float) -> float:
    return CF * (
        -2 * (1 + z) * np.log((1 - z) / z) - 4 * np.log(z) / (1 - z) + 6 + 4 * z
    )


def cq_NLO_S(z: float) -> float:
    return CF * (4.0 * np.log(1.0 - z) - 3.0) / (1.0 - z)


def cq_NLO_L(x: float) -> float:
    return CF * (
        4.0 * np.log(1.0 - x) ** 2.0 / 2.0
        - 3.0 * np.log(1.0 - x)
        - (9.0 + 4.0 * zeta(2))
    )


cq_NLO = RSL(regular=cq_NLO_R, singular=cq_NLO_S, local=cq_NLO_L)


# NLO gluon
def cg_NLO_R(z: float) -> float:
    return (
        (2.0 - 4.0 * z * (1.0 - z)) * np.log((1.0 - z) / z) - 2.0 + 16.0 * z * (1.0 - z)
    )


def cg_NLO_S(z: float) -> float:
    return 0.0


def cg_NLO_L(x: float) -> float:
    return 0.0


cg_NLO = RSL(regular=cg_NLO_R, singular=cg_NLO_S, local=cg_NLO_L)


# Okay, now we have all necessary formulae and codes - let's do it! Lets implement a simple DIS code using the DIS factorization formula and the coefficient functions from above.
#
# Actually, for the sake of this tutorial, let's *not* do this in one big step, but in three smaller steps. By the way, this is a good advice for programming: if you can break it in smaller steps, you should break it in smaller steps.

# ### Excercise 1

# Let's abstract the task of convoltion an RSL coefficient function with another function first.
#
# Write a function to convolute an `RSL` object with an other function.

# In[4]:


def convolute(c: RSL, f: Callable[[float], float], x: float) -> float:
    r"""Convolute c(z) with f(z) with respect to x: :math:`(c\otimes f)(x)`."""
    # write your code here ...
    # recall the RSL definition from above!
    # at some point you want to use quad for the integration ...
    # see https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.quad.html#scipy.integrate.quad
    return 0.0


# In[6]:


# Check your implementation here!
np.testing.assert_allclose(
    convolute(
        RSL(regular=lambda _z: 1.0, singular=lambda _z: 0.0, local=lambda _x: 0.0),
        lambda z: z * (1.0 - z),
        0.5,
    ),
    0.125,
)
np.testing.assert_allclose(
    convolute(
        RSL(regular=lambda z: z, singular=lambda _z: 0.0, local=lambda _x: 0.0),
        lambda z: z * (1.0 - z),
        0.5,
    ),
    0.09657359027972266,
)
np.testing.assert_allclose(
    convolute(
        RSL(regular=lambda _z: 0.0, singular=lambda _z: 0.0, local=lambda _x: 1.0),
        lambda z: z * (1.0 - z),
        0.5,
    ),
    0.25,
)
np.testing.assert_allclose(
    convolute(
        RSL(
            regular=lambda _z: 0.0,
            singular=lambda z: 1.0 / (1.0 - z),
            local=lambda _x: 0,
        ),
        lambda z: z * (1.0 - z),
        0.5,
    ),
    0.04828679513973628,
)
np.testing.assert_allclose(
    convolute(
        RSL(
            regular=lambda z: z,
            singular=lambda z: 1.0 / (1.0 - z),
            local=lambda x: np.log(1.0 - x),
        ),
        lambda z: z * (1.0 - z),
        0.5,
    ),
    -0.0284264097205274,
)


# In[7]:


# We have a small mathematical hickup in our NLO quark coefficient function:
# a removable singularity at z = 1 through log(z)/(1-z) .
# For details see also https://en.wikipedia.org/wiki/Removable_singularity
np.testing.assert_allclose(
    convolute(
        RSL(
            regular=lambda z: np.log(z) / (1.0 - z),
            singular=lambda _z: 0.0,
            local=lambda _x: 0.0,
        ),
        lambda z: z * (1.0 - z),
        0.5,
    ),
    -0.1431167583560283,
)


# ### Excercise 2

# Next, let's take a look at the flavor structure: in our master formula we wrote
# $$ \frac 1 x F_2(x,Q^2) = \sum_{a} C_a \otimes f_a^p $$
# where the sum $a$ runs over all the 5 light quarks ($d,\bar d, u, \bar u, s,\bar s, c, \bar c, b, \bar b$) and the gluon $g$ (so 11 elements).
# However, at NLO there actually only 2 different type of input configuration possible: a quark or a gluon.
# The information which quark is entering is only needed for the coupling (i.e. the respective electric charge) which is a mere prefactor
# and in a similar way we can factorize for the gluon channel the information on how many quarks are coupling to the photon.
# This means to NLO we can write
# $$ \frac 1 x F_2(x,Q^2) = \sum_{b=\{q,g\}} c_b \otimes l_b^p $$
# where we refer to the functions $l_b$ as luminosity channel. The precise mapping between $C_a$ and $c_a$ is also given above.
#
# Define the two luminosity functions $l_q^p(x,Q^2)$ and $l_g^p(x,Q^2)$. (For the scope of this tutorial, we are not computing *grids* here, but use directly a given PDF.)

# In[8]:


# Let's use NNPDF4.0 - so we load it via LHAPDF
pdf = lhapdf.mkPDF("NNPDF40_nnlo_as_01180", 0)
# recall that you can access PDFs via
xgluon_at_x0p5_Q2100GeV2 = pdf.xfxQ2(21, 0.5, 100.0)  # = x*g(x=0.5,Q^2=100 GeV^2)
# where the function signature of xfxQ2 is pid,x,Q^2
# and the pid of the gluon is 21 while it is 1=d,-1=dbar,2=u, etc. see also https://pdg.lbl.gov/2022/reviews/contents_sports.html
# also recall that xfxQ2 returns x * PDF, but in our factorization formula we need the true PDF


# recall the deifitions of C and c from above ...
def lumi_q(x: float, Q2: float) -> float:
    """Quark luminosity at momentum fraction x and factorization scale Q2."""
    # write your code here ...
    return 0.0


def lumi_g(x: float, Q2: float) -> float:
    """Gluon luminosity at momentum fraction x and factorization scale Q2."""
    # write your code here ...
    return 0.0


# In[10]:


# Check your implementation here!
np.testing.assert_allclose(lumi_q(0.1, 10.0), 4.294867071979338)
np.testing.assert_allclose(lumi_g(0.1, 10.0), 5.0 * 3.2411771900340334)


# Knowing and identifying the luminosity structure of an observable is one of the major concerns of a theoretical physicist.
# We are always trying to find new observables which expose a certain luminosity structure to access different physics related to different flavors.

# ### Excercise 3

# Finally, we can join everything together.
#
# Write a function to compute a DIS structure function.

# In[11]:


# recall that you can access the strong coupling via
alpha_s_at_100GeV2 = pdf.alphasQ2(100.0)


def f2(x: float, Q2: float) -> float:
    """DIS F2 structure function at Bjorken x and virtuality Q2."""
    # write your code here ...
    # make use of 1 and 2
    # remember to mask the Q2 dependence of 2 for 1 - people refer to this as "Currying" https://en.wikipedia.org/wiki/Currying
    # also recall the normalization of the lhs of our factorization formula
    return 0.0


# In[13]:


# Check you implementation here!
# if you want to start with a leading order comparison (which is a good idea!), uncomment the following line
# np.testing.assert_allclose(f2(0.001, 10), 0.9948344436236175, rtol=3e-4)
# else, we want of course the most precise calculation! (we will compare to real data in the next subsection)
# (Again, let's not worry about IntegrationWarnings for the moment ...)
with warnings.catch_warnings():
    warnings.filterwarnings("ignore", category=IntegrationWarning)
    np.testing.assert_allclose(f2(0.001, 10), 1.0189, rtol=3e-4)
    np.testing.assert_allclose(f2(0.001, 1000), 2.5299, rtol=3e-4)
    np.testing.assert_allclose(f2(0.1, 100), 0.4023, rtol=3e-4)


# ## Compare to data

# Now, we are ready to confront these theoretical calculation with some real life experimental data:
# DIS was and is a major key stone in determining PDFs. Some of the most influential datasets are coming from
# [HERA](https://en.wikipedia.org/wiki/HERA_(particle_accelerator)) which was operating at [DESY](https://en.wikipedia.org/wiki/DESY) in Germany between 1992 and 2007 (there is some data still published nowadays!).
# Here, we use the combined measurement of the two experiments [H1](https://en.wikipedia.org/wiki/H1_(particle_detector)) and [ZEUS](https://en.wikipedia.org/wiki/ZEUS_(particle_detector)) for the collision of an electron with an proton at a center of mass energy of $\sqrt{s} = 318\,\text{GeV}$ available from [hepdata](https://doi.org/10.17182/hepdata.68951.v1/t5) or the [HERA website](https://www.desy.de/h1zeus/herapdf20/) (in practice, here, we use the latter).
# Download the data from
#   https://www.desy.de/h1zeus/herapdf20/HERA1+2_NCem.dat
# for example on UNIX you can run
# $ curl https://www.desy.de/h1zeus/herapdf20/HERA1+2_NCem.dat -o HERA1+2_NCem.dat
#
# This data is also used inside the NNPDF PDF determination, so we should better find some agreement.

# In[14]:


# Load the experimental data via pandas
data = pd.read_csv("./HERA1+2_NCem.dat", sep=r"\s+")
# the (for us) relevant columns are Q2, x and Sigma, which is the actual measurement
pd.concat([data["Q2"], data["x"], data["Sigma"]], axis=1)


# In[15]:


# let's compute our own predictions for the given points!
with warnings.catch_warnings():
    warnings.filterwarnings("ignore", category=IntegrationWarning)
    f2_data = data.apply(
        lambda dat: pd.Series(
            # we can keep the name 'Sigma' for the experimental value and just add a column 'F2' for our numbers
            [dat["x"], dat["Q2"], dat["Sigma"], f2(dat["x"], dat["Q2"])],
            index=["x", "Q2", "Sigma", "F2"],
        ),
        axis=1,
    )


# In[16]:


# let's plot the experimental data against our predictions!
plt.close()
fig, axes = plt.subplots(6, 4, figsize=(10, 15), clear=True, layout="tight")
fig.suptitle(r"HERA $e^-\! + p \,\to\, e^-\! + X$ with $\sqrt{s}=318 \mathrm{GeV}$")
# the data is organised by Q2 - so let's make a plot for each value
for q2, ax in zip(f2_data["Q2"].unique(), axes.flatten()):
    # select the correct data
    dat = f2_data[f2_data["Q2"] == q2]
    ax.set_title(f"$Q^2 = {q2}\\,\\mathrm{{GeV}}^2$")
    # and plot
    (pexp,) = ax.plot(dat["x"], dat["Sigma"], "o", label="exp")
    (pth,) = ax.plot(dat["x"], dat["F2"], "x", label="th")
    # it is convenient to plot log(x) instead of just x
    ax.set_xscale("log")
fig.legend(handles=[pexp, pth], ncols=2)
fig.savefig("HERA1+2_NCem.pdf")


# Using a simple NLO calculation we can already see that we can get some decent agreement between theory and experiment which can give us some confidence that something is working.
# However, we also see that we do not match everywhere - so, instead of saying "a simple NLO calculation" we should say "a naive NLO calculation", because the full story behind DIS is much more complicated - see lecture!
