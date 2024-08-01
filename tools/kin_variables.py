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
def kin_ArgusPsR(p_tracks_sig, p_tracks_tag):
    '''
    Calculates Argus normalized variable
    
    :param p_b: tag side momentum
    :param p_a: signal side momentum
    :param Ecms: CMS energy
    :param M_mom: mom mass
    :return: Argus normalized variable
    '''
    #accepted = True
    #CMS_E = Ecms
    CMS_E = 10.58
    E_mom = CMS_E/2.0
    M_mom = 1.77709 # GeV, Belle II reported mass for the Tau Lepton 

    pa = TLorentzVector()
    pb = TLorentzVector()
    pa_boost = TLorentzVector()
    boostPS = TVector3()

    pa = p_tracks_sig
    pb = p_tracks_tag

    P_mom = np.sqrt(E_mom*E_mom - M_mom*M_mom)
    betav = P_mom/E_mom

    # Boosting of pion in the signal side to the Tau restframe with following approximation
    # a) the momentum direction of the  in the 1-prong side as the opposite direction of the momentum of the lepton in the three-prong side
    # b) the B energy by ECMS/2
    boostPS = betav*(pb.BoostVector().Unit())
    pa.Boost(boostPS)

    Ea_boost = pa.E()

    # Argus normalized variable
    Xval = 2.0 * Ea_boost / M_mom


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

def cos_theta_th(p_tracks_tag):
    '''
    Calculates cos(theta)_tau_hadron
    
    :param p_tracks_tag: tag side momentum

    :return: cos_theta_th
    '''
    
    E_mom = 10.58
    M_mom = 1.77709
    # Tag energy (hadrons energy in the tag side)
    E_tag = p_tracks_tag.E()
    # Tag momentum (hadrons momentum in the tag side)
    p_tag = p_tracks_tag.P()

    p_mom = np.sqrt(E_mom*E_mom - M_mom*M_mom)

    p_nu_tag = E_mom - E_tag

    cos_th = (p_mom**2 + p_tag**2 - p_nu_tag**2) / (2.0 * p_mom * p_tag) 

    return cos_th

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

def add_ArgusPsR(data):

    Xval_list = []

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

        #Ecms = float(data.at[row, "Ecms"])
        #Xval = kin_ArgusPsR(p_tracks_sig=pa, p_tracks_tag=pb, Ecms=Ecms)
        Xval = kin_ArgusPsR(p_tracks_sig=pa, p_tracks_tag=pb)
        Xval_list.append(Xval)
        
    data["Xps"] = Xval_list
    return data

def add_cos_theta_th(data):
    cos_th_list = []
    for row in range(len(data.index)):
        
        p_tracks_tag = TLorentzVector()

        p_tracks_tag.SetPxPyPzE(float(data.at[row, "track_tag_px_CMS"]),
                                float(data.at[row, "track_tag_py_CMS"]),
                                float(data.at[row, "track_tag_pz_CMS"]),
                                float(data.at[row, "track_tag_E_CMS"]))

        cos_th = cos_theta_th(p_tracks_tag)
        cos_th_list.append(cos_th)

    data["cos_th"] = cos_th_list

    return data