class Constants:

    #------------------------constants for monocrystal------------------------------------#
    nbdom = 6
    AS = 1e-5

    # Déformation ferro-électrique à saturation
    epz0=2e-3
    #C/m2 % Polarisation initiale dans un domaine (à champ électrique macro nul)
    Pol0=0.3
    # Elastic coefficients
    C11=106e9
    C12=62e9
    C44=44e9
    # Piezoelectric coefficients (3=polarisation direction)
    D31=-2.1e-10
    D33=4.5e-10
    D15=5.8e-10
    # Dielectric permittivitty (3=polarisation direction)
    k33=2.00e-8
    k11=2.00e-8


    #------------------------constants for polycrystal------------------------------------#
    #grain orientations
    n1 = 7
    n2 = 4
    n3 = 3

    #------------------------constants for script_RC------------------------------------#   
    #for saving plot images
    saveDir="plots/"
    dpi=200
    lineWidth = 1

    #piezo coefficient
    alpha = 1e-3;
    nbt = 16
    nbe = 16
    tSpace = -150
    eSpace = 5