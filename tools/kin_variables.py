import ROOT
from ROOT import TH1D, TVector3, TLorentzVector, TGraph, TMath, TFeldmanCousins
import numpy as np
import pandas as pd
import random

# Mmin2 and Mmax2 calculations from kinematics
def kin_Mmin_Mmax(mom_sig,mom_tag,Ecms):
    '''
    Calculates Mmin2 and Mmax2
    
    :param mom_sig: signal side momentum
    :param mom_tag: tag side momentum
    :param Ecms: CMS energy
    :return: Mmin2 and Mmax2
    '''
    
    accepted = True
    CMS_E = Ecms
    tau_m = 1.77686
    
    x0 = (tau_m/CMS_E)*(tau_m/CMS_E)
    
    pa = TLorentzVector()
    pb = TLorentzVector()

    pa = mom_sig
    pb = mom_tag
    n_a = TVector3()
    n_b = TVector3()
    n_a = (1.0/CMS_E)*pa.Vect()
    n_b = (1.0/CMS_E)*pb.Vect()
    ab = n_a.Dot(n_b)
        
    za = pa.E()/CMS_E
    zb = pb.E()/CMS_E
        
    a2 = n_a.Mag2()
    b2 = n_b.Mag2()
        
    wa = (zb*zb - zb - b2 - 2*ab)
    wb = (za*za - za + a2)
        
    H = TVector3()
    H = wa*n_a + wb*n_b; 

    acrossb = TVector3()
    acrossb = n_a.Cross(n_b)
            
    A1 = b2
    A2 = a2
    A3 = 2.0*ab
    B1 = 2.0*(n_b.Dot(H))
    B2 = 2.0*(n_a.Dot(H))
    C1 = 4.0*acrossb.Mag2()
    D1 = H.Dot(H) -C1*(0.5 - za)*(0.5 - za)
    #A = A1 + A2 + A3
    #B = B1 + B2
    #C = D1
    A = A1
    B = -B1 + C1 -2*A1*x0 - A3*x0
    C = A3*x0*x0 + A2*x0*x0 + A1*x0*x0 + B2*x0 + B1*x0 + D1
    
    
    if (B*B - 4.0*A*C) < 0 : 
        accepted = False
        Mmin2 = 1000
        Mmax2 = 1000
    else :
        accepted = True
        Mmax2 = CMS_E*CMS_E*(-B + np.sqrt(B*B - 4*A*C))/(2*A)
        Mmin2 = CMS_E*CMS_E*(-B - np.sqrt(B*B - 4*A*C))/(2*A)  
        
    return accepted, Mmin2, Mmax2

# Smearing in momenta
def P_smearing(P, perc):
    '''
    Applies smearing in momenta
    
    :param P: momentum
    :param perc: percentage of smearing
    :return: smeared momentum
    '''

    # Let's create the detector effect (smearing 5% around in momenta)
    pX_CMS = P.Px()
    pY_CMS = P.Py()
    pZ_CMS = P.Pz()
    p_ECMS = P.Energy()

    pt_cms = (pX_CMS**2 + pY_CMS**2 + pZ_CMS**2)**0.5 # Calculate transverse momentum
    sigma_pt = perc * pt_cms  # Standard deviation of 1% of transverse momentum

            # Applying smearing using a Gaussian distribution with mean = 1 and standard deviation = sigma_pt
    pX_CMS_smeared = pX_CMS * (1.0 + random.gauss(0,sigma_pt))
    pY_CMS_smeared = pY_CMS * (1.0 + random.gauss(0,sigma_pt))
    pZ_CMS_smeared = pZ_CMS * (1.0 + random.gauss(0,sigma_pt))

    return pX_CMS_smeared, pY_CMS_smeared, pZ_CMS_smeared, p_ECMS


# Argus normalized variable calculation
def kin_ArgusPsR(p_b,p_a, Ecms, M_mom):
    '''
    Calculates Argus normalized variable
    
    :param p_b: tag side momentum
    :param p_a: signal side momentum
    :param Ecms: CMS energy
    :param M_mom: mom mass
    :return: Argus normalized variable
    '''
    #accepted = True
    CMS_E = Ecms
    E_B = CMS_E/2.0
    M_mom = M_mom

    pa = TLorentzVector()
    pb = TLorentzVector()
    pa_boost = TLorentzVector()
    boostPS = TVector3()

    pa = p_a
    pb = p_b

    P_B = np.sqrt(E_B*E_B - M_mom*M_mom)
    betav = P_B/E_B

    # Boosting of pion in the signal side to the B restframe with following approximation
    # a) the momentum direction of the B in the 1-prong side as the opposite direction of the momentum of the D0+l in the three-prong side
    # b) the B energy by ECMS/2
    boostPS = betav*(pb.BoostVector().Unit())
    pa.Boost(boostPS)

    Ea_boost = pa.E()

    Xval = (2.0*Ea_boost)/M_mom

    return Xval

# Q2 calculation
def kin_q2(Ecms, px, py, pz):
    """
    Analisis based on B+ -> K+ nunu
    Formula found in next link, page: 28:
    https://docs.belle2.org/record/3785/files/BELLE2-TALK-DRAFT-2023-117.pdf
    :param Ecms: CMS energy
    :param px: px
    :param py: py
    :param pz: pz
    :return: q2
    """

    # Mass of the K
    m_sig = 0.493677
    s = Ecms**2
    
    pt_cms = (px**2 + py**2 + pz**2)**0.5
    Esig_CMS = (m_sig**2 + pt_cms**2)**0.5

    q2 = s/4 + m_sig**2 - np.sqrt(s)*Esig_CMS

    return q2

# Adding Mmin2 and Mmax2 to the dataframe
def add_Mmin2Mmax2(data):
    '''
    Adds Mmin2 and Mmax2 to the dataframe
    those variables could be found in the paper:
    https://journals.aps.org/prd/abstract/10.1103/PhysRevD.102.115001
    
    :param data: dataframe
    :return: dataframe with Mmin2 and Mmax2 added
    '''

    Mmin2_list = []
    Mmax2_list = []

    for row in range(len(data.index)):

        pa = TLorentzVector()
        pb = TLorentzVector()

        pa.SetPxPyPzE(float(data.at[row, "track_sig_px_CMS"]),
                      float(data.at[row, "track_sig_py_CMS"]),
                      float(data.at[row, "track_sig_pz_CMS"]),
                      float(data.at[row, "track_sig_E_CMS"]))

        pb.SetPxPyPzE(float(data.at[row, "track_tag_px_CMS"]),
                      float(data.at[row, "track_tag_py_CMS"]),
                      float(data.at[row, "track_tag_pz_CMS"]),
                      float(data.at[row, "track_tag_E_CMS"]))

        Ecms = float(data.at[row, "Ecms"])
        flag, Mmin2, Mmax2 = kin_Mmin_Mmax(pa,pb,Ecms)
        Mmin2_list.append(float(Mmin2))
        Mmax2_list.append(float(Mmax2))
    
    # Using 'Mmin2' and 'Mmax2' as the column name and equating it to the list
    data["Mmin2"] = Mmin2_list
    data["Mmax2"] = Mmax2_list

    return data 