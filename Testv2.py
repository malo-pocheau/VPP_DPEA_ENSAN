# -*- coding: utf-8 -*-
"""
Created on Mon May  9 14:37:32 2022

@author: Cyprien Courtois
"""
import math as math
import numpy as np
import matplotlib.pyplot as plt

#%%
''' Appendices '''
''' 
ENTREE : 
    Paramètres extérieurs :
        Derive (angle de dérive)
        rho eau
        Vb (vitesse navire) en kts
        Tcq (Tirant d'eau coque)
        T (Tirant d'eau avec quille)
        Heeling (angle de gîte)
        GM 
        Visco (viscosité de l'eau)
    Paramètres de géometrie de Quille:
        Vk (déplacement)
        Cordepied (Corde point le plus bas)
        Cordebout (corde niveau coque)
        Xleadquille (position X bord d'attaque pt plus haut)
        Yleadquille (position Y bord d'attaque pt plus haut)
        Zleadquille (position Z bord d'attaque pt plus haut)
        H (hauteur quille)
        Tc (rapport t/c quille, t epaisseur)
        Cplaque (Coefficient d'effet de plaque de quille)
        Calagegite (angle de calage en gite de quille)
        Calageincidence (calage incidence quille)
        Alpha0 (coefficient alpha 0 des profils asymétriuqes de quille6)
        Lg_bulbe (longueur de bulbe)
        Tc_bulbe (t/c bulbe)
        
    Paramètres de géometrie de Safran :
        nb_saf (nombre de safran)
        cordepiedsaf (corde safran pt plus bas)
        cordeboutsaf (corde safran niveau coque)
        Xaxesaf (Postiion X bord d'attaque safran pt plus haut)
        Yaxesaf
        Zaxesaf
        Hsaf (hauteur safran)
        Tcsaf (rapport t/c safran, t epaisseur)
        Cplaquesaf (coefficient d'effet de plaque safran)
        Calagegitesaf (angle de calage en gite safran)
        Calageincidencesaf (calage d'incidence safran)
        Alphasaf (coefficient alpha 0 des profils asymétriques safran)
        
    
SORTIE :
    Fx
    Fy
    Zce
    Mx
    Pour :
        Quille
        Bulbe
        Safran 1
        Safran 2
        Total
'''
#%%
''' 
Transformation d'angle de degres au radians
'''

def degtorad(a):
    return(a*math.pi/180)

'''
Transformation d'angle de radians en degres
'''

def radtodeg(a):
    return(a*180/math.pi)

#%%

''' Second test pour Appendices
Creation d'une fonction unique pour appendices avec différenciation interne
'''

# Type : Quille, Bulbe, MonoSafran, BiSafran
# Centrage : Centre, Auvent, Sousvent (centre par défault, utile pour bisafran)

def Appendice (Type="None", Centrage="Centre", Cordepied=0, Cordebout=0, Calagegite=0, Tc=0, Tcb=0, H=0, Calageincidence=0, Alpha0=0, Heeling=0,  Vb=0, Cplaque=0, Visco=0, rho=0, Derive=0, Tcq=0, T=0, Zlead=0, L=0, YAxe=0, ZAxe=0, GM=0):
    if Type == "None":
        return("Entrer un type d'appendice")
        
    if Type == "Quille":
        k = 1
    else :
        k = 0.9
    
    Moycorde = (Cordepied + Cordebout)/2
    kforme = 2*Tcq+60*Tcq**4
    Re = Vb/1.944*Moycorde/Visco
    Pdyn = 0.5*rho*(k*Vb/1.944)**2
    Sproj = Moycorde*H
    Sw = 2.2*Sproj
    Angle_profilsym = Calageincidence - Alpha0
    AoA = Derive + Angle_profilsym
    
    if Type == "Bulbe" : 
        Re = Vb/1.944*L/Visco
        Cdi = (1+kforme)*0.075/((math.log10(Re)-2)**2)
        Sw = 2*math.pi*L*(0.5*Tcb/L)        
        Drag = Cdi*Pdyn*Sw
        Fx = Drag*math.cos(degtorad(AoA))
        Fy = Drag*math.sin(degtorad(AoA))
        Zce = (Zlead - (H+L*Tcb))*math.cos(degtorad(Heeling))
        Mx = Fy*Zce/math.cos(degtorad(Heeling))
        
        Data = [Re, kforme, Cdi, Sw, Pdyn, Drag]
        Name = ["Re = ", "kforme = ", "Cdi = ", "Sw = ", "Pdyn = ", "Drag = "]
        Result = ["Bulbe"]
        for i in range (0,len(Data)) :
            Result += [Name[i]+str(Data[i])]
        return(Fx, Fy, Zce, Mx, Result)
    
    
    
    Hproj_heel = H*math.cos(degtorad(Heeling + Calageincidence))
    Sproj_heel = Hproj_heel*Moycorde
    Apc = Cplaque * (Hproj_heel**2)/Sproj_heel
    
    Sweep = radtodeg(math.atan((Cordebout - Cordepied)/H))
        
    Cl = 5.7*2*Apc/(1.8 + math.cos(degtorad(Sweep))*(((2*Apc)**2)/math.cos(degtorad(Sweep))**4 + 4)**0.5)*degtorad(Derive)
    Lift = 0.5*rho*(k*Vb/1.944)**2*Sproj*Cl
    C_hull = 1.8*(Tc/(-Zlead+H))+1
    C_heel = 1 - 0.382*degtorad(Heeling)
    Lift_eff = Lift *C_hull * C_heel
    Cdi = Cl**2/(math.pi*4*Apc)
    Cdv = (1 + kforme)*0.075/((math.log10(Re)-2)**2)
    Di = Cdi*Pdyn*Sproj_heel
    Dv = Cdv*Pdyn*Sw
    Drag = Di + Dv
    Cgeom = (H/3)*(Cordebout+2*Cordepied)/(Cordepied + Cordebout)
    
    
    Fx = Lift_eff*math.sin(degtorad(AoA))*math.cos(degtorad(Heeling))-Drag*math.cos(degtorad(AoA))
    Fy = Lift_eff*math.cos(degtorad(AoA))*math.cos(degtorad(Heeling))+Drag*math.sin(degtorad(AoA))
    
        

    
    if Type == "Quille" : 
        Zce = (Zlead-Cgeom)*math.cos(radtodeg(Heeling))
    
    if Type == "MonoSafran":
        Zce = (ZAxe-Cgeom)*math.cos(radtodeg(Heeling))
    
    if Type == "BiSafran" :
        if Centrage == "Auvent" :
            gite = Calagegite + Heeling
            Angle_boutsaf = radtodeg(math.atan(YAxe/(-ZAxe+GM)))
            Angle_piedsaf = radtodeg(math.atan((YAxe + H*math.sin(degtorad(Calagegite)))/(GM-ZAxe+H)))
            CGgoem = 0.5*H*Cordepied+(H*(Cordebout-Cordepied)/6)
            if Heeling > (90-Angle_piedsaf):    
                Liftcorr = 0
                Dragcorr = 0
            else :
                if Heeling < (90-Angle_boutsaf):
                    Liftcorr = Lift_eff
                    Dragcorr = Drag
                else :
                    Liftcorr = Lift_eff*(1-(Heeling-Angle_piedsaf)/(Angle_boutsaf-Angle_piedsaf))
                    Dragcorr = Drag*(1-(Heeling-Angle_piedsaf)/(Angle_boutsaf-Angle_piedsaf))
            Fx = Liftcorr*math.sin(degtorad(AoA))*math.cos(degtorad(Heeling))-Dragcorr*math.cos(degtorad(AoA))
            Fy = Liftcorr*math.cos(degtorad(AoA))*math.cos(degtorad(Heeling))+Dragcorr*math.sin(degtorad(AoA))
            if Heeling < (90 - Angle_boutsaf):
                
                Zce = (ZAxe-Cgeom)*math.cos(degtorad(AoA))
            else :
                Zce = 0
        if Centrage == "Sousvent" :
            gite = -Calagegite + Heeling
            Zce = (Zlead-Cgeom)*math.cos(radtodeg(gite))
    Mx = Fy*Zce/math.cos(degtorad(Heeling))
    
    Data = [Moycorde, kforme, Sproj, Sw, Angle_profilsym, Hproj_heel, Sproj_heel, Apc, Re, Pdyn, AoA, Sweep, Cl, Lift, C_hull, C_heel, Lift_eff, Cdi, Cdv, Di, Dv, Drag, Cgeom]
    Name = ["Moycorde = ", "kforme = ", "Sproj = ", "Sw = ", "Angle_profilsym = ", "Hproj_heel = ", "Sproj_heel = ", "Apc = ", "Re = ", "Pdyn = ", "AoA = ", "Sweep = ", "Cl = ", "Lift = ", "C_hull = ", "C_heel = ", "Lift_eff = ", "Cdi = ", "Cdv = ", "Di = ", "Dv = ", "Drag = ", "Cgeom = "]
    
    Result = []
    for i in range (0,len(Data)) :
        Result += [Name[i]+str(Data[i])]
    
    return(Fx, Fy, Zce, Mx, Result)

#Tests

print("Quille", Appendice(Type ="Quille", Cordepied=1, Cordebout=1.5, Tc=0.4, H=1.9, Calageincidence=0, Alpha0=0, Heeling=0.1,  Vb=10, Cplaque=2, Visco=10**(-6), rho=1025, Derive=0.1, Tcq=0.15, T=2.3, Zlead=-0.4))
print("Bulbe", Appendice(Type="Bulbe", Vb=10, L=2, Visco=10**(-6), Tcb=0.4, rho=1025, Zlead=-0.4, H=1.9, Heeling=0.1))
print("MonoSafran", Appendice (Type="MonoSafran", Cordepied=0.3, Cordebout=0.5, Vb=10, YAxe=1, ZAxe=-0.2, H=0.5, Tc=0.15, Cplaque=1.6, Calagegite=30, Calageincidence=0, Visco=10**(-6), rho=1025, Heeling=0.1, Derive=0.1))
print("BiSafran au vent", Appendice(Type="BiSafran", Centrage = "Auvent", Cordepied=0.3, Cordebout=0.5, Vb=10, YAxe=1, ZAxe=-0.2, H=0.5, Tc=0.15, Cplaque=1.6, Calagegite=30, Calageincidence=0, Visco=10**(-6), rho=1025, Heeling=0.1, Derive=0.1))
print("BiSafran Sous vent", Appendice(Type="BiSafran", Centrage = "Sousvent", Cordepied=0.3, Cordebout=0.5, Vb=10, YAxe=1, ZAxe=-0.2, H=0.5, Tc=0.15, Cplaque=1.6, Calagegite=30, Calageincidence=0, Visco=10**(-6), rho=1025, Heeling=0.1, Derive=0.1))


#%%
'''
Fonction avant fusion 

'''

'''
def quille (Cordepied, Cordebout, Tc, H, Calageincidence, Alpha0, Heeling,  Vb, Cplaque, Visco, rho, Derive, Tcq, T, Zleadquille):
    Moycorde = (Cordepied + Cordebout)/2
    kforme = 2*Tcq+60*Tcq**4
    Sproj = Moycorde*H
    Sw = 2.2*Sproj
    Angle_profilsym = Calageincidence - Alpha0
    Hproj_heel = H*math.cos(degtorad(Heeling + Calageincidence))
    Sproj_heel = Hproj_heel*Moycorde
    Apc = Cplaque * (Hproj_heel**2)/Sproj_heel
    Re = Vb/1.944*Moycorde/Visco
    Pdyn = 0.5*rho*(Vb/1.944)**2
    AoA = Derive + Angle_profilsym
    Sweep = radtodeg(math.atan((Cordebout - Cordepied)/H))
        
    Cl = 5.7*2*Apc/(1.8 + math.cos(degtorad(Sweep))*(((2*Apc)**2)/math.cos(degtorad(Sweep))**4 + 4)**0.5)*degtorad(Derive)
    Lift = 0.5*rho*(Vb/1.944)**2*Sproj*Cl
    C_hull = 1.8*(Tc/(-Zleadquille+H))+1
    C_heel = 1 - 0.382*degtorad(Heeling)
    Lift_eff = Lift *C_hull * C_heel
    Cdi = Cl**2/(math.pi*4*Apc)
    Cdv = (1 + kforme)*0.075/((math.log10(Re)-2)**2)
    Di = Cdi*Pdyn*Sproj_heel
    Dv = Cdv*Pdyn*Sw
    Drag = Di + Dv
    Cgeom = (H/3)*(Cordebout+2*Cordepied)/(Cordepied + Cordebout)
    
    Data = [Moycorde, kforme, Sproj, Sw, Angle_profilsym, Hproj_heel, Sproj_heel, Apc, Re, Pdyn, AoA, Sweep, Cl, Lift, C_hull, C_heel, Lift_eff, Cdi, Cdv, Di, Dv, Drag, Cgeom]
    Name = ["Moycorde = ", "kforme = ", "Sproj = ", "Sw = ", "Angle_profilsym = ", "Hproj_heel = ", "Sproj_heel = ", "Apc = ", "Re = ", "Pdyn = ", "AoA = ", "Sweep = ", "Cl = ", "Lift = ", "C_hull = ", "C_heel = ", "Lift_eff = ", "Cdi = ", "Cdv = ", "Di = ", "Dv = ", "Drag = ", "Cgeom = "]
    Fx = Lift_eff*math.sin(degtorad(AoA))*math.cos(degtorad(Heeling))-Drag*math.cos(degtorad(AoA))
    Fy = Lift_eff*math.cos(degtorad(AoA))*math.cos(degtorad(Heeling))+Drag*math.sin(degtorad(AoA))
    Zce = (Zleadquille-Cgeom)*math.cos(radtodeg(Heeling))
    
    Result = []
    for i in range (0,len(Data)) :
        Result += [Name[i]+str(Data[i])]
    
    return(Fx, Fy, Zce, Result)


Fx, Fy, Zce, Result = quille (Cordepied=1, Cordebout=1.5, Tc=0.4, H=1.9, Calageincidence=0, Alpha0=0, Heeling=0.1,  Vb=10, Cplaque=2, Visco=10**(-6), rho=1025, Derive=0.1, Tcq=0.15, T=2.3, Zleadquille=-0.4)
print("Cas de quille")
print(Result)
print("Fx = ", Fx, "Fy = ", Fy, "Zce = ", Zce)

#Correct vis a vis de l'excel sauf : Zce --> Rad ou Deg ? 


def bulbe (Vb, L, visco, Tcb, rho, AoA, Zleadquille, H, Heeling):
    Re = Vb/1.944*L/visco
    kforme = 2*Tcb+60*Tcb**4
    Cdi = (1+kforme)*0.075/((math.log10(Re)-2)**2)
    Sw = 2*math.pi*L*(0.5*Tcb/L)
    Pdyn = 0.5*rho*(Vb/1.944)**2
    
    Drag = Cdi*Pdyn*Sw
    
    Data = [Re, kforme, Cdi, Sw, Pdyn, Drag]
    Name = ["Re = ", "kforme = ", "Cdi = ", "Sw = ", "Pdyn = ", "Drag = "]
    
    Fx = Drag*math.cos(degtorad(AoA))
    Fy = Drag*math.sin(degtorad(AoA))
    Zce = (Zleadquille - (H+L*Tcb))*math.cos(degtorad(Heeling))
            
    Result = []
    for i in range (0,len(Data)) :
        Result += [Name[i]+str(Data[i])]
        
    return(Fx, Fy, Zce, Result)

Fx, Fy, Zce, Result = bulbe (Vb=10, L=2, visco=10**(-6), Tcb=0.4, rho=1025, AoA=0.1, Zleadquille=-0.4, H=1.9, Heeling=0.1)
print("Cas de Bulbe")
print(Result)
print("Fx = ", Fx, "Fy = ", Fy, "Zce = ", Zce)
'''

 
'''

def monosafran():
    Moycorde =
    kforme =
    Sproj =
    Sw =
    Hproj_heel =
    Sproj_heel =
    Apc =
    Re =
    Aoa =
    Sweep =
    
    Moycorde = (Cordepied + Cordebout)/2
    kforme = 2*Tcq+60*Tcq**4
    Sproj = Moycorde*H
    Sw = 2.2*Sproj
    Angle_profilsym = Calageincidence - Alpha0
    Hproj_heel = H*math.cos(degtorad(Heeling + Calageincidence))
    Sproj_heel = Hproj_heel*Moycorde
    Apc = Cplaque * (Hproj_heel**2)/Sproj_heel
    Re = Vb/1.944*Moycorde/Visco
    Pdyn = 0.5*rho*(Vb/1.944)**2
    AoA = Derive + Angle_profilsym
    Sweep = radtodeg(math.atan((Cordebout - Cordepied)/H)) 
    
    Cl =
    Lift =
    Lift_eff =
    Cdi =
    Cdv =
    Di =
    Dv =
    Drag =
    Cgeom =
    
    Fx =
    Fy =
    Zce =
    return()

def bisafran():
    Moycorde =
    Gite =
    kforme =
    Sproj =
    Sw =
    Hproj_heel =
    Sproj_heel =
    Apc =
    Re =
    Aoa =
    Sweep =

    Cl =
    Lift =
    Lift_eff =
    Cdi =
    Cdv =
    Di =
    Dv =
    Drag =
    Cgeom =
    
    if n = 0:
        Fx =
        Fy =
        Zce =
    
    if n = 1 :
        Angle_boutsaf =
        Angle_piedsaf =
        Cg_geom =
        Lift_corr =
        Drag_corr =
        Fx =
        Fy = 
        Zce =
        

    return()

'''
    
