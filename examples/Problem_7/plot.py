import h5py
import matplotlib.pyplot as plt
import numpy as np

f = h5py.File('output.h5', 'r');

phi = f['scalar_flux'];
z = f['z'];

ref1 = np.array([0.1184152951,0.1528865341,0.1746178112,0.1899757963,0.2015841429,0.210674174,0.2179227253,0.2237596943,0.2284874497,0.2323318467,0.2354671823,0.2380303724,0.2401300901,0.2418531801,0.2432693924,0.2444349745,0.2453954423,0.2461877392,0.2468419311,0.2473825456,0.2478296351,0.2481996241,0.2485059889,0.2487598032,0.2489701787,0.2491446217,0.2492893224,0.2494093903,0.2495090465,0.2495917806,0.249660479,0.2497175311,0.2497649151,0.2498042696,0.2498369522,0.249864088,0.2498866091,0.2499052876,0.249920763,0.249933564,0.2499441277,0.2499528139,0.2499599184,0.2499656829,0.2499703031,0.2499739359,0.2499767045,0.2499787029,0.2499799987,0.2499806358])
ref2 = np.array([0.295127282,0.3899724402,0.4620655386,0.5214566878,0.5724664161,0.6171774357,0.6567607783,0.6919869198,0.723428437,0.7515446227,0.7767199744,0.7992839537,0.8195226296,0.8376864423,0.8539958713,0.8686458216,0.8818091419,0.8936395126,0.9042738545,0.9138343622,0.922430237,0.9301591772,0.9371086673,0.9433571004,0.9489747615,0.9540246901,0.9585634429,0.9626417671,0.9663051987,0.9695945931,0.9725465987,0.9751940772,0.9775664797,0.9796901816,0.9815887814,0.9832833673,0.9847927545,0.9861336967,0.987321074,0.9883680595,0.9892862664,0.990085878,0.9907757609,0.9913635643,0.9918558047,0.9922579389,0.9925744245,0.9928087694,0.9929635712,0.9930405455])
ref3 = np.array([        0.448932795,0.5969456841,0.7153719629,0.8172906326,0.9082368055,0.9907939194,1.0663538559,1.1357991516,1.1997730714,1.2587907222,1.3132878497,1.3636444725,1.4101978289,1.4532504138,1.4930755124,1.5299212985,1.5640140206,1.5955605613,1.6247505393,1.6517580683,1.6767432508,1.699853463,1.7212244723,1.7409814191,1.7592396866,1.7761056764,1.7916775046,1.8060456297,1.8192934206,1.8314976721,1.8427290747,1.8530526421,1.8625281011,1.8712102478,1.8791492725,1.8863910559,1.8929774387,1.898946467,1.9043326149,1.9091669855,1.9134774923,1.9172890216,1.920623577,1.9235004073,1.9259361182,1.9279447683,1.929537951,1.9307248615,1.93151235,1.9319049613])
ref4 = np.array([        1.0574651844,1.4164973152,1.7223378609,2.0001270124,2.2601337479,2.5068085524,2.7422788498,2.9677103139,3.1838455074,3.3912231177,3.5902715407,3.7813518465,3.9647795125,4.1408366627,4.3097796698,4.4718442462,4.6272490324,4.77619821,4.9188834495,5.0554853832,5.1861747373,5.311113211,5.4304541708,5.544343205,5.6529185739,5.7563115809,5.8546468831,5.9480427557,6.0366113192,6.1204587383,6.1996853976,6.2743860583,6.3446499994,6.4105611453,6.4721981819,6.5296346631,6.5829391075,6.6321750879,6.677401312,6.7186716968,6.7560354366,6.7895370633,6.8192165024,6.8451091211,6.8672457726,6.885652833,6.9003522335,6.9113614876,6.9186937116,6.9223576411])

plt.plot(z,phi, 'o')
plt.plot(z,ref1)
plt.plot(z,ref2)
plt.plot(z,ref3)
plt.plot(z,ref4)
plt.xlabel("z, cm");
plt.ylabel("scalar flux, /cm^2.s");
plt.show(); 
