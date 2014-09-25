import numpy as np
import math

def SNPHWE(obs_hets, obs_hom1, obs_hom2):
    """
    Calculate exact Hardy-Weinberg equilibrium.
    Directly adapted from R-code at: http://www.sph.umich.edu/csg/abecasis/Exact/r_instruct.html
    Citation: 
        Wigginton JE, Cutler DJ and Abecasis GR
        Am J Hum Genet (2005) 76: 887-93

    Input:
                  HET   HOM1  HOM2
        MARKER_1  100    3    5
        MARKER_2   57   14   50
        MARKER_3   31   32   51
        MARKER_4   47    3    5

    Output:
        MARKER_1   6.19509447581e-21
        MARKER_2   0.842279756571
        MARKER_3   2.56377346576e-06
        MARKER_4   1.16298486151e-07
    """
    if (obs_hets < 0) or (obs_hom1 < 0) or (obs_hom2 < 0):
        return -1
    
    # total number of genotypes
    N = obs_hets + obs_hom1 + obs_hom2

    # rare homozygotes, common homozygotes
    obs_homr = min(obs_hom1, obs_hom2)
    obs_homc = max(obs_hom1, obs_hom2)

    # number of rare allele copies
    rare = (2 * obs_homr) + obs_hets

    # Initialize probability array
    probs = np.zeros(rare + 1, dtype=np.float)

    mid = math.floor(rare * ( 2 * N - rare) / (2 * N))
    if (mid % 2) != (rare % 2):
        mid += 1
    
    probs[mid] = 1.0
    mysum = 1.0

    curr_hets = mid
    curr_homr = (rare - mid) / 2
    curr_homc = N - curr_hets - curr_homr
    while curr_hets >= 2:
        probs[curr_hets - 2] = probs[curr_hets] * curr_hets * (curr_hets - 1.0) / (4.0 * (curr_homr + 1.0)  * (curr_homc + 1.0))
        mysum += probs[curr_hets - 2]
        # 2 fewer heterozygotes -> add 1 rare homozygote, 1 common homozygote
        curr_hets -= 2
        curr_homr += 1
        curr_homc += 1

    curr_hets = mid
    curr_homr = (rare - mid) / 2
    curr_homc = N - curr_hets - curr_homr
   
    while curr_hets <= (rare - 2):
        probs[curr_hets + 2] = probs[curr_hets] * 4.0 * curr_homr * curr_homc / ((curr_hets + 2.0) * (curr_hets + 1.0))
        mysum += probs[curr_hets + 2]
        # add 2 heterozygotes -> subtract 1 rare homozygtote, 1 common homozygote
        curr_hets += 2
        curr_homr -= 1
        curr_homc -= 1
    target = probs[obs_hets]
    return min(1.0, np.sum(probs[probs <= target]) / float(mysum))