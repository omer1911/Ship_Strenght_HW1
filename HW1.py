
import numpy as np
from scipy.interpolate import interp1d
from scipy.integrate import cumtrapz
import matplotlib.pyplot as plt
offset = np.array([
        [0.0, 0.0, 0.0, 0.0, 0.0, 0.09, 0.324],
        [0.011, 0.055, 0.071, 0.083, 0.114, 0.335, 0.604],
        [0.047, 0.131, 0.173, 0.218, 0.314, 0.581, 0.813],
        [0.203, 0.343, 0.458, 0.573, 0.772, 0.98, 1.094],
        [0.439, 0.632, 0.809, 0.978, 1.089, 1.151, 1.187],
        [0.693, 0.908, 1.085, 1.18, 1.2, 1.2, 1.199],
        [0.856, 1.049, 1.183, 1.2, 1.2, 1.2, 1.2],
        [0.768, 0.993, 1.151, 1.188, 1.193, 1.191, 1.193],
        [0.502, 0.788, 0.986, 1.066, 1.087, 1.101, 1.123],
        [0.211, 0.484, 0.672, 0.751, 0.787, 0.82, 0.863],
        [0.041, 0.194, 0.295, 0.341, 0.357, 0.374, 0.429],
        [0.013, 0.084, 0.135, 0.149, 0.155, 0.165, 0.207],
        [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.024]])*10     # hidro yarı genişlik değerleri 
boy = 120
genislik = 24
draft = 8.73 
yogunluk = 1.025
D=13.09
Cb=0.61
def bonjeanAlani(offset, posta0, posta, suhatti):
    offset_new = np.zeros((201, 7))
    for i in range(7):
        f = interp1d(posta0, offset[:, i], kind = "cubic")
        offset_new[:, i] = f(posta)
    alan = np.zeros((201, 7))   # BON-JEAN ALANLARI
    for i in range(201):
        alan[i, 1:] = 2 * cumtrapz(offset_new[i, :], suhatti[:])
    return offset_new,alan
def prohaskaDagilimi(boy, deplasman):
    qx = np.zeros(201)   # TOPLAM GEMİ AĞIRLIK DAĞILIMI (PROHASKA YÖNTEMİ)
    x=(deplasman*-1.8672453217741563/(boy**2))*54/7 #ilk durumda ki sakin su LCG-LCB arası fark değeri girilmiş.
    a = (.68 * deplasman / boy)+x
    b = 1.185 * deplasman / boy
    c = (.58 * deplasman / boy)-x
    n = 40
    posta=np.linspace(0,120,201)
    for i in range(201):
        if i < 68:
            qx[i] = a + posta[i] * (b - a) / 40   #uzunluğa bölündü oran-orantı yapıldı
        elif i > 135:
            qx[i] = b - n * (b - c) / 40  #uzunluğa bölündü oran-orantı yapıldı
            n -= 0.6
    qx[67:136]=b   
    qx[136:]=qx[136:][::-1]  #listeye n ilk 40 değeri ile girdiği için c 136. eleman olarak geldi bu yüzden liste ters çevrildi.
    return qx
def ataletDagilimi(boy, posta):
    Ix = np.zeros(201)   # ATALET MOMENT DAĞILIMI
    Wmin =4.032 #minimum midship section modulus [GL LOYD SECTİON 5C 2.1]
    Iy = 3 * Wmin * boy / 100  # ORTA KESİT ATALET MOMENTİ [m4]
    for i in range(201):
        if posta[i] <= boy / 20:
            Ix[i] = 5 * Iy * posta[i] / boy
        if boy / 20 < posta[i] <= 7 * boy / 20:
            Ix[i] = .25 * Iy + (15 * Iy) * (posta[i] - boy / 20) / (6 * boy)
        if 7 * boy / 20 < posta[i] <= 15 * boy / 20:
            Ix[i] = Iy
        if 15 * boy / 20 < posta[i] <= 19 * boy / 20:
            Ix[i] = Iy - 2.5 * (Iy / boy) * (posta[i] - 15 * boy / 20)
        if posta[i] > 19 * boy / 20:
            Ix[i] = .5 * Iy - 10 * (Iy / boy) * (posta[i] - 19 * boy / 20)
    return Ix
def hesap_su(name,ax,qx,posta):
    px = ax-qx
    dpx0 = np.array([0, *cumtrapz(px,posta)])
    '''print(dpx0)
    print(max(abs(dpx0))*0.03)'''
    dpx = dpx0 - dpx0[-1] * posta / boy# LİNEER DÜZELTME (son değer en büyük değerin 0.03 unden küçük)
    ddpx0 = np.array([0, *cumtrapz(dpx,posta)]) # LİNEER DÜZELTME (son değer en büyük değerin 0.06 sından küçük)
    '''print(ddpx0)
    print(max(abs(ddpx0))*0.06)'''
    if name == "cukur" or   name =='sakin' or name=='tepe':
            ddpx = ddpx0 - ddpx0[-1] * posta / boy 
            Qx = dpx*boy/200
            Mx = ddpx*(boy/200)**2

    return Qx,Mx
def hesaplamalar(offset, boy, draft, yogunluk):
    # Define your helper functions here if they are not already implemented
    def grafik(name, ax, qx, Qx, Mx,Stress):
        plt.figure(figsize=(10, 4))
        plt.title(name, fontweight="bold")
        plt.xlabel("Gemi Boyu [m]", fontweight="bold")
        plt.plot(posta, ax, posta, qx, posta, Qx, posta, Mx,posta,Stress)
        plt.grid(color="red", linestyle="-", linewidth=.3)
        plt.legend(["a [ton/m]", "q [ton/m]", "Q/5 [ton]", "M/60 [ton.m]","Stress[Mpa]"], loc="best")
        plt.show()
    # Define your other helper functions (e.g., bonjeanAlani, ataletDagilimi, alan_eq, prohaskaDagilimi, hesap_su, dalgaKaydirma) here
    # Constants and arrays
    posta0 = np.array([0, .5, 1, 2, 3, 4, 5, 6, 7, 8, 9, 9.5, 10]) * boy / 10
    posta = np.linspace(0, boy, 201)
    suhatti = np.array([0, .3, 1, 2, 3, 4, 5]) * draft / 4
    offset_new, alan = bonjeanAlani(offset, posta0, posta, suhatti)
    Ix = ataletDagilimi(boy, posta)
    polinomlar = alan_eq(suhatti, alan)
    # Calculation for calm water
    ax = alan[:, 5] * yogunluk
    deplasman_new = np.trapz(ax, posta)
    qx = prohaskaDagilimi(boy, deplasman_new)
    qx_n = np.trapz(qx, posta)
    LCG = np.sum(qx * posta) / np.sum(qx)
    LCB = np.sum(alan[:, 5] * posta) / np.sum(alan[:, 5])

    # Print results for calm water
    print("Real Cb =", Cb)
    print("Calc. Disp =", deplasman_new)
    print("Calc. Cb =", deplasman_new / (boy * draft * genislik))
    print("LCB(SAKİN) =", LCB)
    print("LCG(SAKİN) =", LCG)
    ax_n = np.trapz(ax, posta)
    print("Deplasman(Ax)(SAKİN) =", ax_n)
    print("Deplasman(qx)(SAKİN) =", qx_n)
    # Create and show the graphs for calm water
    Qx, Mx = hesap_su("sakin", ax, qx, posta)
    # Additional calculations for calm water
    max_m = max(abs(Mx))
    max_q = max(abs(Qx))
    max_m_index = np.argmax(np.abs(Mx))
    max_q_index = np.argmax(np.abs(Qx))
    mx_m_from_q = np.abs(Mx[max_q_index])
    # Print additional results for calm water
    print("Maksimum_Momentyeri_Sakin_su =", max_m_index)
    print("Maksimum_kesmeyeri_Sakin_su =", max_q_index)
    print("Maksimum_Moment_Sakin_su =", max_m)
    print("Maksimum_kesme_Sakin_su =", max_q)
    print("Max Kesmedeki Moment =", mx_m_from_q)
    E_gerilmesi=max_m*(D-7.2)/Ix[max_m_index]*9.81/1000    #Maksimum eğilme gerilmesi kayma gerilmesi 0 kabul edildi #momentin maximum olduğu posta 102 gelmiştir. 102. postanın alan merkezi aşağıdaki kordinatlar(suhattina denk gelen alan değerlerim) ile bulundu 
                                                           #(0,0),(12.5139,0.65475),(46.69,2.1825),(98.73,4.365),(151.106,10.9125),(203.47,8.73),(255.85,10.9125),        Cross Sectional Properties
                                                           #(-255.85,10.9125),(-203.47,8.73),(-151.106,10.9125),(-98.73,4.365),(-46.69,2.1825),(-12.5139,0.65475),(0,0)    Ömer Can Demir Github(https://github.com/omer1911/cross_section_properties)
                                                           #sakin su için momentin maximum olduğu postanın alanının N.A'sı=7.2 metre gelmiştir
    print("Maksimum Eğilme Gerilmesi Sakin Su","=",E_gerilmesi) 
                                                          #kesmenin maximum olduğu posta 158 gelmiştir. 158. postanın alan merkezi aşağıdaki kordinatlar(suhattina denk gelen alan değerlerim) ile bulundu
                                                          #(0,0),(4.92,0.65475),(23.62931876,2.1825),(56.31,4.365),(91.55,10.9125),(128.47,8.73),(166.2,10.9125),      
                                                           #(-166.2,10.9125),(-128.47,8.73),(-91.55,10.9125),(-56.31,4.365),(-23.62931876,2.1825),(-4.92,0.65475),(0,0)                                
    E_gerilmesi_t=mx_m_from_q*(D-7.37-1)/Ix[max_q_index]*9.81/1000 ##sakin su için kesmenin maximum olduğu postanın alanının N.A'sı=7.37 metre gelmiştir. ymax=D-7.37-1m aşağısı alınacak
    S=polinomlar[max_q_index](D)-polinomlar[max_q_index](D-1)*5.22 #daha onceden dalga kaydırma için oluşrutulan alan değerlerinin denklemi burda D ve D nin 1 m aşağıdaki alanı çıkarılarak 1.ALAN MOMENT değeri nulundu
    kayma_gerilmesi=(max_q*S/(Ix[max_q_index]*np.sqrt(boy)))*9.81/1000
    Von_mises=np.sqrt(E_gerilmesi_t**2+3*kayma_gerilmesi**2)
    print("Sakin Su Von mises","=",Von_mises)
    all_stress=np.zeros(201)
    for i in range(201):
        all_stress[i]=Mx[i]*(D-7.26)/Ix[i]*9.81/1000                                                        
    grafik("sakin", ax, qx, Qx / 5, Mx / 60,all_stress) 
    
    # Calculation for wavy water in a cavity
    ax, Alan = dalgaKaydirma("cukur", boy, draft, deplasman_new, alan, posta, suhatti)
    LCB = np.sum(Alan * posta) / np.sum(Alan)
    Qx, Mx = hesap_su("cukur", ax, qx, posta)
    # Print results for wavy water in a cavity
    print("LCB(CUKUR) =", LCB)
    print("LCG(CUKUR) =", LCG)
    ax_n = np.trapz(ax, posta)
    print("Deplasman(Ax)(CUKUR) =", ax_n)
    print("Deplasman(qx)(CUKUR) =", qx_n)
    # Create and show the graph for wavy water in a cavity
    
    # Additional calculations for wavy water in a cavity
    max_m = max(abs(Mx))
    max_q = max(abs(Qx))
    max_m_index = np.argmax(np.abs(Mx))
    max_q_index = np.argmax(np.abs(Qx))
    mx_m_from_q = np.abs(Mx[max_q_index])
    # Print additional results for wavy water in a cavity
    print("Maksimum_Momentyeri_cukur =", max_m_index)
    print("Maksimum_kesmeyeri_cukur =", max_q_index)
    print("Maksimum_Moment_cukur =", max_m)
    print("Maksimum_kesme_cukur =", max_q)
    print("Max Kesmedeki Moment =", mx_m_from_q)
    E_gerilmesi=max_m*(D-7.26)/Ix[max_m_index]*9.81/1000                #  momentin maximum olduğu posta 130 gelmiştir. 130. postanın alan merkezi aşağıdaki kordinatlar(suhattina denk gelen alan değerlerim) ile bulundu 
    all_stress=np.zeros(201)
    original_array = np.linspace(0, 120, 201)

# Y değerleri belirleniyor (ilk ve son değer sabit, ortadaki değer 7.26)
    y_values = np.array([0, 7.26, 0])

# Interpolasyon yapılıyor
    interpolated_array = np.interp(posta, [original_array[0], original_array[-1]], y_values)
    for i in range(201):
        all_stress[i]=Mx[i]*(D-7.26)/Ix[i]*9.81/1000                                                        #(0,0),(10.17,0.65475),(40.63,2.1825),(89.44,4.365),(139.85,10.9125),(190.58,8.73),(241.61,10.9125),      
    grafik("cukur", ax, qx, Qx / 5, Mx / 60,all_stress)                                                                 #(-241.61,10.9125),(-190.58,8.73),(-139.85,10.9125),(-89.45,4.365),(-40.63,2.1825),(-10.17,0.65475),(0,0) 
                                                                       # ymax=7.26m geldi.Q kesme küçük olduğu için ihmal edildi.
    print("Maksimum Eğilme Gerilmesi cukur ","=",E_gerilmesi) 
    E_gerilmesi_t=mx_m_from_q*(D-7.43-1)/Ix[max_q_index]*9.81/1000 ##cukur için kesmenin maximum olduğu postanın(173. posta) alanının N.A'sı=7.43 metre gelmiştir. ymax=D-7.43-1m aşağısı alınacak
    S=polinomlar[max_q_index](D)-polinomlar[max_q_index](D-1)*5.16 #daha onceden dalga kaydırma için oluşrutulan alan değerlerinin denklemi burda D ve D nin 1 m aşağıdaki alanı çıkarılarak 1.ALAN MOMENT değeri bulundu
    kayma_gerilmesi=(max_q*S/(Ix[max_q_index]*np.sqrt(boy)))*9.81/1000
    Von_mises=np.sqrt(E_gerilmesi_t**2+3*kayma_gerilmesi**2)
    print("Dalgali-Cukur Von mises","=",Von_mises)
    # Similar calculations for wavy water at the peak
    ax, Alan = dalgaKaydirma("tepe", boy, draft, deplasman_new, alan, posta, suhatti)
    LCB = np.sum(Alan * posta) / np.sum(Alan)
    Qx, Mx = hesap_su("tepe", ax, qx, posta)
    # Print results for wavy water at the peak
    print("LCB(TEPE) =", LCB)
    print("LCG(TEPE) =", LCG)
    ax_n = np.trapz(ax, posta)
    print("Deplasman(Ax)(TEPE) =", ax_n)
    print("Deplasman(qx)(TEPE) =", qx_n)
    # Create and show the graph for wavy water at the peak
    # Additional calculations for wavy water at the peak
    max_m = max(abs(Mx))
    max_q = max(abs(Qx))
    max_m_index = np.argmax(np.abs(Mx))
    max_q_index = np.argmax(np.abs(Qx))
    mx_m_from_q = np.abs(Mx[max_q_index])
    E_gerilmesi=max_m*(D-7.2)/Ix[max_m_index]*9.81/1000 
    print("Maksimum Eğilme Gerilmesi Tepe","=",E_gerilmesi)      # maximum moment 101. postada geldi ve 102. posta ile paralel posta olarak kabul yapıldı. dolayısıyla sakin su için  hesaplanan NA=7.2 max eğilme gerilmesi için aynen geçerli.
                                  
                                   
                                                                 #(0,0),(6.10,0.65475),(25.62931876,2.1825),(60.75,4.365),(102.6,10.9125),(149.2,8.73),(199.85,10.9125),      
                                                                 #(-199.85,10.9125),(-149.2,8.73),(-102.6,10.9125),(-60.75,4.365),(-25.62931876,2.1825),(-6.10,0.65475),(0,0)
    E_gerilmesi_t=mx_m_from_q*(D-7.45-1)/Ix[max_q_index]*9.81/1000 ##tepe durmur için kesmenin maximum olduğu postanın alanının N.A'sı=7.45 metre gelmiştir. ymax=D-7.45-1m aşağısı alınacak
    S=polinomlar[max_q_index](D)-polinomlar[max_q_index](D-1)*5.14 #daha onceden dalga kaydırma için oluşrutulan alan değerlerinin denklemi burda D ve D nin 1 m aşağıdaki alanı çıkarılarak 1.ALAN MOMENT değeri bulundu
    kayma_gerilmesi=(max_q*S/(Ix[max_q_index]*np.sqrt(boy)))*9.81/1000
    Von_mises=np.sqrt(E_gerilmesi_t**2+3*kayma_gerilmesi**2)
    print("Dalgalı-Tepe Von mises Tepe","=",Von_mises)  
    all_stress=np.zeros(201)
    for i in range(201):
        all_stress[i]=Mx[i]*(D-7.26)/Ix[i]*9.81/1000                                                        
    grafik("tepe", ax, qx, Qx / 5, Mx / 60,all_stress) 
def dalgaKaydirma(name,boy, draft, deplasman,  alan, posta, suhatti ):
    H=boy/20
    yogunluk=1.025
    polinomlar=alan_eq(suhatti,alan)
    if name =="cukur":
        dalga_katsayi = [1, .966, .871, .795, .578, .422, .28, .16, .072, .018, 0]
        posta1 = np.array([0, .5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5]) * boy / 10
        c = np.interp(posta[:101], posta1, dalga_katsayi)   # YARI DALGA DEĞERLERİ
        c = np.concatenate((np.delete(c, -1), np.flipud(c)), axis = 0) # YARI DALGAYI TAM YAPMA
        ksi  = draft - H / 2 + H * c
        ksi = np.roll(ksi, -5) #dalgayı yatay kaydırma
    elif name=="tepe":
        dalga_katsayi = [0, .018, .072, .16, .28, .422, .578, .795, .871, .966,1]
        posta1 = np.array([0, .5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5]) * boy / 10
        c = np.interp(posta[:101], posta1, dalga_katsayi)   # YARI DALGA DEĞERLERİ
        c = np.concatenate((np.delete(c, -1), np.flipud(c)), axis = 0) # YARI DALGAYI TAM YAPMA
        ksi = draft - H / 2 + H * c
        ksi = np.roll(ksi, -1)
    ax = np.zeros(201)   
    alan=np.zeros(201)
    for i in range(200):
        for j in range (201):
            alan[j]=polinomlar[j](ksi[j])
            ax[j] = yogunluk * polinomlar[j](ksi[j])
        if round(deplasman) < round(np.trapz(ax, posta)):  # DALGAYI DİKEY KAYDIRARAK a(x) BULMAK
            ksi -= 0.0077
        elif round(deplasman) > round(np.trapz(ax, posta)):
            ksi += 0.0077
        else: break 
    return ax,alan

def alan_eq(suhatti,alan):
    polinomlar = np.zeros(201, dtype=object)  # Her bir postadaki alanlar için denklem oluşturulması
    degree = 6  # 6 tane su hatti ile yapıldığı için derece 6 olarak belrilendi
    for i in range(201):
        coefficients = np.polyfit(suhatti, alan[i], degree)
        polinomlar[i] = np.poly1d(coefficients)
    return polinomlar
hesaplamalar(offset,boy,draft,yogunluk)


    
