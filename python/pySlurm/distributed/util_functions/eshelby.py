import numpy as np

def Eshelby33_muphylin(C_inf, form_inclu):
    if form_inclu == 1:
        Ndemag33=1/3*np.eye(3)
    elif form_inclu == 2:
        Ndemag33=np.matrix([[1/2, 0, 0],
                            [0, 1/2, 0],
                            [0, 0, 0]])

    else:
        raise Exception("Eshelby 33 form_inclu value not valid")

    return Ndemag33

def Eshelby66_muphylin(C_inf, icalc):
    if icalc == 1:

        nu0=C_inf.item(0,1)/(C_inf.item(0,0)+C_inf.item(0,1));

        JJ=(1/3)*np.concatenate((np.concatenate((np.ones((3,3)), np.zeros((3,3))), axis=1), \
            np.concatenate((np.zeros((3,3)), np.zeros((3,3))), axis=1)), axis=0)

        KK=(1/3)*np.concatenate((np.concatenate((-1*np.ones((3,3))+3*np.eye(3,3), np.zeros((3,3))), axis=1), \
            np.concatenate((np.zeros((3,3)), 3*np.eye(3,3)), axis=1)), axis=0)

        SEsh66=(1+nu0)/(3*(1-nu0))*JJ+(2/15)*(4-5*nu0)/(1-nu0)*KK;

    elif icalc == 2:
        raise Exception("Eshelby 66 icalc case not coded")

    elif icalc == 3:
        raise Exception("Eshelby 66 icalc case not coded")

    elif icalc == 4:
        raise Exception("Eshelby 66 icalc case not coded")

    else:
        raise Exception("Eshelby 66 icalc value not valid")

    return SEsh66