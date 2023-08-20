# -*- coding: UTF-8 -*-
import numpy as np
import math as mt
from datetime import datetime
from tkinter import *
from tkinter import filedialog
import os
from tkinter import messagebox

def Title():
    title = [
#str(datetime.now())+'\n\n',
'Transitivity - Code\n\n',
#'Program developed by:\nHugo Gontijo Machado, Flávio Olimpio Sanches Neto,\nNayara Dantas Coutinho, Kleber Carlos Mundim,\nVincenzo Aquilanti and Valter Henrique Carvalho Silva.\n\n',
]
    return title

class Extract:
    pi = float((mt.atan2(1.0, 1.0)) * 4)
    h = 6.62607004E-34#J s
    kb = 1.3806503E-23 #J/K
    T = [250.0, 273.15, 298.15, 300.0, 310.0]
    Tinv1000 = [float(1000.0/x) for x in T]
    Tinv = [float(1.0/x) for x in T]
    LnT = [float(np.log(x)) for x in T]
    LogT = [float(np.log10(x)) for x in T]
    def __init__(self,filename,type,DoFcp):
        self.p = {}

        self.p['DoFcp'] = DoFcp
        self.p['filename'] = filename
        self.p['zpe'] = 0.0 ; self.p['En'] = 0.0 ; self.p['Ent'] = 0.0 ; self.p['EnG'] = 0.0
        self.p['DoF'] = 0 ; self.p['freqi'] = 0.0 ; self.p['Tvib'] = [] ; self.p['Trot'] = [] ;self.p['NSym'] = 0.0
        self.p['mass'] = 0.0 ; self.p['masskg'] = 0.0 ; self.p['P'] =0.0

        self.p['type'] = type
        self.p['Qe'] = 0.0 ; self.p['Qt'] = [] ; self.p['Qr'] = [] ; self.p['Qv'] = [] ; self.p['QTot'] = []

        self.p['rVol'] = 0.0

        self.p['title'] = Title()

    def ext(self):
        file = open(self.p['filename'],'r', encoding="utf-8")
        for ln in file:
            if "Deg. of freedom" in ln:
                if self.p['DoFcp'] == -1:
                    self.p['DoF'] = int(ln.split()[3])
                elif self.p['DoFcp'] >= 0:
                    self.p['DoF'] = self.p['DoFcp']
            if "Kelvin.  Pressure" in ln:
                self.p['P'] = float(ln.split()[4]) * 101325
            if "Molecular mass:" in ln:
                self.p['mass'] = float(ln.split()[2])
                self.p['masskg'] = self.p['mass'] * 1.660538921E-27
            if "Rotational symmetry number" in ln or 'ROTATIONAL SYMMETRY NUMBER' in ln:
                self.p['NSym'] = float(ln.split()[3])
            if "Rotational temperature" in ln:
                self.p['Trot'].extend([float(tr) for tr in ln.split()[3:]])
            if "Vibrational temperatures:" in ln:
                self.p['Tvib'].extend([float(tv) for tv in ln.split()[2:]])
                line = file.readline()
                self.p['Tvib'].extend([float(tv) for tv in line.split()[1:]])
                line = file.readline()
                while line.strip() != '':
                    self.p['Tvib'].extend([float(tv) for tv in line.split()])
                    line = file.readline()
            if "Zero-point correction" in ln:
                self.p['zpe'] = (float(ln.split()[2])) * 627.509
                __zpe = (float(ln.split()[2]))
            if "Sum of electronic and zero-point Energies" in ln:
                self.p['En'] = float(ln.split()[6]) * 627.509
            if "Sum of electronic and thermal Enthalpies" in ln:
                self.p['Ent'] = (float(ln.split()[6]) + __zpe) * 627.509
            if "Sum of electronic and thermal Free Energies" in ln:
                self.p['EnG'] = (float(ln.split()[7]) + __zpe) * 627.509
            if "Electronic      0" in ln:
                self.p['Qe'] = float(ln.split()[1].replace('D','E'))
            if self.p['type'] == "ts":
                if "frequencies ---" in ln:
                    try: errorts += 1
                    except: self.p['freqi'] = -float(ln[20:30]) ; errorts = 1
            if "a0 for SCRF calculation" in ln:
                self.p['rVol'] =  float( ln.split()[6] )
        file.close()

    def Atom(self):

        if self.p['mass'] == 0.0 or self.p['Qe'] == 0.0 or self.p['En'] == 0.0:
            file = open(self.p['filename'], 'r', encoding="utf-8")
            iso=0
            for ln in file:
                if 'iso=' in ln:
                    iso = int(ln.split()[0].split('=')[1].split(')')[0])
                if 'Input orientation:' in ln:
                    for i in range(5):
                        line = file.readline()
                    atom_n = int(line.split()[1])
                if 'Multiplicity =' in ln:           #For atoms, the electronic partition coefficient is the multiplicity
                    ln = ln.split()
                    self.p['Qe'] = float(ln[ln.index('Multiplicity')+2])
                if 'SCF Done:' in ln:
                    self.p['En'] = float(ln.split()[4]) * 627.509
            file.close()

            dicc = {
                1:1.00783,2:4.00260,3:7.01600,4:9.01218,5:11.00931,6:12.00000,7:14.00307,8:15.99491,9:18.99840,10:19.99244,
                11:22.98977,12:23.98505,13:26.98154,14:27.97693,15:30.97376,16:31.97207,17:34.96885,18:39.96238,
                19:38.96371,20:39.96259,21:44.95591,22:47.94795,23:50.94396,24:51.94051,25:54.93805,26:55.93494,27:58.93320,
                28:57.93535,29:62.92960,30:63.92915,31:68.92558,32:73.92118,33:74.92160,34:79.91652,35:78.91834,36:83.91151,
                37:84.91170,38:87.90560,39:88.90540,40:89.90430,41:92.90600,42:97.90550,43:98.90630,44:101.90370,45:102.90480,46:105.90320,47:106.90509,48:113.90360,49:114.90410,50:117.90180,51:120.90380,52:129.90670,53:126.90040,54:131.90420
            }
            iso_dicc = {
                1:{2:2.01410,3:3.01605},
                6:{13:13.00335,14:14.00324},
                7:{15:15.00011},
                8:{17:16.99913,18:17.99916},
                17:{37:36.96590}
            }

            if iso == 0:
                self.p['mass'] = dicc[atom_n]
            else:
                try:
                    self.p['mass'] = iso_dicc[atom_n][iso]

                except:
                    self.p['mass'] = float(iso)
            self.p['masskg'] = self.p['mass'] * 1.660538921E-27
        if self.p['P'] == 0.0:
            self.p['P'] = 1.0 * 101325E0
        if self.p['Ent'] == 0 or self.p['EnG'] == 0:
        #Enthalpy
            R = (1.9872E-3)  # kcal/(mol*K)
            T = 298.15  # K
            kb = 0.0019872041  # kcal/(mol*K)
            Hcorr = (((3 / 2) * R * T) + (kb * T))  # kcal/mol
            self.p['Ent'] = self.p['En'] + Hcorr # kcal/mol

        #En. Gibbs
            qt = ( ( (2.0 * self.pi * self.p['masskg'] * self.kb * T) / (self.h**2.0) ) ** 1.50 ) * ( ( self.kb * T ) / self.p['P'] )
            qe = self.p['Qe']
            St = R * (np.log(qt) + 1 + (3 / 2))
            Se = R * (np.log(qe))
            Stot = St + Se # kcal/(mol*K)
            Gcorr = Hcorr - (T * Stot)  # kcal/mol

            self.p['EnG'] = self.p['En'] + Gcorr # kcal/mol
    def calc_Q(self):
       #----------------------------- Calculo Q Vibracional-----------------------------------------------------------------
        if self.p['DoF'] == 0:
            for j in range(len(self.T)):
                self.p['Qv'].append(1.0)
        else:
            for j in range(len(self.T)):
                vibT = [-Tvib for Tvib in self.p['Tvib']]
                self.p['Qv'].append(1.0 / (np.product(1 - np.exp(np.divide(vibT, self.T[j])))))
    #----------------------------- Calculo Q Translacional---------------------------------------------------------------
        for j in range(len(self.T)):
            self.p['Qt'].append( ( ( (2.0 * self.pi * self.p['masskg'] * self.kb * self.T[j]) / (self.h**2.0) ) ** 1.50 ) * ( ( self.kb * self.T[j] ) / self.p['P'] ) )
    # ----------------------------- Calculo Q Rotacional------------------------------------------------------------------
        if self.p['DoF'] == 0:
            for j in range(len(self.T)):
                self.p['Qr'].append(1.0)
        elif self.p['DoF'] == 1:
            for j in range(len(self.T)):
                self.p['Qr'].append( self.T[j] / (self.p['NSym'] * self.p['Trot'][0]) )
        elif self.p['DoF'] > 1:
            for j in range(len(self.T)):
                self.p['Qr'].append( ( (self.pi**0.5) / self.p['NSym'] ) * ( (self.T[j]**1.5) / (np.product(self.p['Trot'][:]))**0.5 ) )
    # ----------------------------- Calculo Qp Total------------------------------------------------------------------
        for j in range(len(self.T)):
            self.p['QTot'].append(self.p['Qt'][j] * self.p['Qr'][j] * self.p['Qv'][j] * self.p['Qe'])


########################################################################################################################
########################################################################################################################

class Rate_Calc:
    h = 6.62607004E-34
    pi = (mt.atan2(1.0, 1.0)) * 4
    hb = h / (2.0 * pi)
    r = 1.9872E0
    kb = 1.3806503E-23
    NA = 6.0221409E23
    m2c = 1.0E6
    c = 2.99793E10
    def __init__(self,reac1,reac2,ts,prod1,prod2):
        self.T = ts.T
        self.Tinv1000 = Extract.Tinv1000
        self.Tinv = Extract.Tinv
        self.LnT = Extract.LnT
        self.LogT = Extract.LogT

        self.reac1 = reac1
        self.reac2 = reac2
        self.ts = ts
        self.prod1 = prod1
        self.prod2 = prod2

        self.p ={}

        self.p['freqi'] =  ts.p['freqi']

        self.p['PES_REAC'] = ((reac1.p['En'] + reac2.p['En']) - (reac1.p['En'] + reac2.p['En']))
        self.p['PES_PROD'] = ((prod1.p['En'] + prod2.p['En']) - (reac1.p['En'] + reac2.p['En']))
        self.p['PES_TS'] = ((ts.p['En']) - (reac1.p['En'] + reac2.p['En']))
        self.p['E0'] = (ts.p['En'] - (reac1.p['En'] + reac2.p['En']))  #kcal
        self.p['E0cal'] = self.p['E0'] *1000 #cal
        self.p['H0'] = (ts.p['Ent'] - (reac1.p['Ent'] + reac2.p['Ent'])) * 1000.0 #Cal
        self.p['de'] = ((prod1.p['En'] + prod2.p['En']) - ( reac1.p['En'] + reac2.p['En']))
        self.p['dh'] = ((prod1.p['Ent'] + prod2.p['Ent']) - ( reac1.p['Ent'] + reac2.p['Ent']))
        self.p['dg'] = ((prod1.p['EnG'] + prod2.p['EnG']) - ( reac1.p['EnG'] + reac2.p['EnG']))
        self.p['Tc'] = (self.c * self.h * ts.p['freqi'] / (2.0 * self.pi * self.kb))

        self.p['d'] = (-0.004769426830E0) * ((self.c * self.h * self.NA * ts.p['freqi']) / self.p['E0cal']) ** 2

        self.p['di'] = 1.0 / self.p['d']
        self.p['alpha'] = (2. * self.pi) / (self.h * self.c * self.p['freqi'])


        self.p['QTot'] = []   ; self.p['ktst'] = []    ;self.p['Lnktst'] = []   ;self.p['kdtst'] = [];self.p['Lnkdtst'] = []
        self.p['bell35'] = [] ; self.p['kbell35'] = [] ;self.p['Lnkbell35'] = []
        self.p['fator'] = []  ; self.p['bell58'] = []  ;self.p['kbell58'] = [] ; self.p['Lnkbell58'] = []
        self.p['bell2T'] = [] ; self.p['kbell2T'] = [] ;self.p['Lnkbell2T'] = []
        self.p['beta'] = []   ; self.p['ST'] = []      ;self.p['kST'] = []     ; self.p['LnkST'] = []
        self.p['rates'] = []  ; self.p['Lnrates'] = [] ;self.p['Logrates'] = []

        self.p['Logktst'] = []     ; self.p['Logkbell58'] = [] ; self.p['Logkdtst'] = [] ; self.p['LogkST'] = []
        self.p['Logkbell35'] = []  ; self.p['Logkbell2T'] = []

        self.p['out'] = [] ; self.p['title'] = []

        self.WriteTitle()
    def Calc(self):
        for j in range(len(self.T)):
            X1 = (self.kb * self.T[j] / self.h)
            X2 = (np.exp(-self.p['E0cal'] / (self.r * self.T[j])))
            #            # --------------------- Cáclulo do Coeficiente de partição Total da reação p/ cada temperatura------------------------
            if self.reac2.p['type'] != 'reactant':
                self.p['QTot'].append(self.ts.p['QTot'][j] / self.reac1.p['QTot'][j])
            else:
                self.p['QTot'].append(self.ts.p['QTot'][j] / (self.reac1.p['QTot'][j] * self.reac2.p['QTot'][j]))
            #            # --------------------- Cáclulo do Taxa convencional (kTST)-----------------------------------------------------------

            if self.reac2.p['type'] != 'reactant':
                self.p['ktst'].append(X1 * (self.p['QTot'][j]) * X2)
            else:
                self.p['ktst'].append(self.m2c * X1 * (self.p['QTot'][j]) * X2 / self.NA )
            try:
                self.p['Lnktst'].append(np.log(self.p['ktst'][j]))
                self.p['Logktst'].append(np.log10(self.p['ktst'][j]))
            except:
                self.p['Lnktst'].append(float("NaN"))
                self.p['Logktst'].append(float("NaN"))
            #            # --------------------- Cáclulo do Taxa d-TST ------------------------------------------------------------------------
            if self.reac2.p['type'] != 'reactant':
                self.p['kdtst'].append(
                    X1 * self.p['QTot'][j] * np.exp(self.p['di'] * np.log(1 - self.p['d'] * self.p['E0cal'] / (self.r * self.T[j]))))
            else:
                self.p['kdtst'].append(self.m2c * X1 * self.p['QTot'][j] * np.exp(self.p['di'] * np.log(1 - self.p['d'] * self.p['E0cal'] / (self.r * self.T[j]))) / self.NA)
            try:
                self.p['Lnkdtst'].append(np.log(self.p['kdtst'][j]))
                self.p['Logkdtst'].append(np.log10(self.p['kdtst'][j]))
            except:
                self.p['Lnkdtst'].append(float("NaN"))
                self.p['Logkdtst'].append(float("NaN"))
            # --------------------- Correção de Bell 1935 -----------------------------------------------------------------------

            self.p['bell35'].append((self.kb * self.T[j] - 1.0E0 * self.hb * self.p['freqi'] * self.c * np.exp(
                -(2.0E0 * 4.180E0) * self.pi * self.p['E0cal'] / (self.c * self.h * self.p['freqi'] * self.NA)) / X2) / (
                                           self.kb * self.T[j] - self.hb * self.p['freqi'] * self.c))
            self.p['kbell35'].append(self.p['bell35'][j] * self.p['ktst'][j])
            try:
                self.p['Lnkbell35'].append(np.log(self.p['kbell35'][j]))
                self.p['Logkbell35'].append(np.log10(self.p['kbell35'][j]))
            except:
                self.p['Lnkbell35'].append(float("NaN"))
                self.p['Logkbell35'].append(float("NaN"))
            # --------------------- Correção de Bell 1958 -----------------------------------------------------------------------

            self.p['fator'].append((0.50E0 * self.h * self.p['freqi'] * self.c) / (self.kb * self.T[j]))
            self.p['bell58'].append((self.p['fator'][j] / np.sin(self.p['fator'][j])))
            self.p['kbell58'].append(self.p['bell58'][j] * self.p['ktst'][j])

            try:
                self.p['Lnkbell58'].append(np.log(self.p['kbell58'][j]))
                self.p['Logkbell58'].append(np.log10(self.p['kbell58'][j]))
            except:
                self.p['Lnkbell58'].append(float("NaN"))
                self.p['Logkbell58'].append(float("NaN"))

            self.p['bell2T'].append(self.p['bell58'][j] - 1.0E0 * np.exp(
                -(2.0E0 * 4.180E0) * self.pi * self.p['E0cal'] / (self.c * self.h * self.p['freqi'] * self.NA)) / (
                                       X2 * (1.0E0 * self.pi / self.p['fator'][j] - 1.0E0)))
            self.p['kbell2T'].append(self.p['bell2T'][j] * self.p['ktst'][j])
            try:
                self.p['Lnkbell2T'].append(np.log(self.p['kbell2T'][j]))
                self.p['Logkbell2T'].append(np.log10(self.p['kbell2T'][j]))
            except:
                self.p['Lnkbell2T'].append(float("NaN"))
                self.p['Logkbell2T'].append(float("NaN"))
            # --------------------- Correção ST ---------------------------------------------------------------------------------

            self.p['beta'].append(1.0E0 / (self.kb * self.T[j]))

            if self.p['beta'][j] <= self.p['alpha']:
                if self.p['dh'] <= 0.0E0:
                    dh = 0.0E0
                    self.p['ST'].append(self.p['bell58'][j] - (self.p['beta'][j] / ((self.p['alpha']) - self.p['beta'][j]) * np.exp(
                        ((self.p['beta'][j] - (self.p['alpha']))) * (self.p['E0cal'] - dh) * (4.1868E0 / self.NA))))
                else:
                    self.p['ST'].append(self.p['bell58'][j] - (self.p['beta'][j] / ((self.p['alpha']) - self.p['beta'][j]) * np.exp(
                        ((self.p['beta'][j] - (self.p['alpha']))) * (self.p['E0cal'] - self.p['dh']) * (4.1868E0 / self.NA))))
            elif self.p['beta'][j] > self.p['alpha']:
                if self.p['dh'] <= 0.0E0:
                    dh = 0.0E0
                    self.p['ST'].append((self.p['beta'][j] / (self.p['beta'][j] - (self.p['alpha']))) * (
                            np.exp(
                                ((self.p['beta'][j] - (self.p['alpha']))) * ((self.p['E0cal'] - dh) * (4.1868E0 / self.NA))) - 1.0))
                else:
                    self.p['ST'].append((self.p['beta'][j] / (self.p['beta'][j] - (self.p['alpha']))) * (
                            np.exp(
                                ((self.p['beta'][j] - (self.p['alpha']))) * ((self.p['E0cal'] - self.p['dh']) * (4.1868E0 / self.NA))) - 1.0))
            self.p['kST'].append(self.p['ST'][j] * self.p['ktst'][j])
            try:
                self.p['LnkST'].append(np.log(self.p['kST'][j]))
                self.p['LogkST'].append(np.log10(self.p['kST'][j]))
            except:
                self.p['LnkST'].append(float("NaN"))
                self.p['LogkST'].append(float("NaN"))
    def WriteTitle(self):
        self.p['title'] = Title()
        if self.reac2.p['type'] != 'reactant':
            self.p['title'].append("# REACTION:  " + (self.reac1.p['filename'].split('/')[-1]).split('.')[0] + "  --->  " +
                         (self.prod1.p['filename'].split('/')[-1]).split('.')[0] + "\n\n")
        else:
            self.p['title'].append("# REACTION:  " + (self.reac1.p['filename'].split('/')[-1]).split('.')[0] + " + " +
                            (self.reac2.p['filename'].split('/')[-1]).split('.')[0] + "  --->  " +
                            (self.prod1.p['filename'].split('/')[-1]).split('.')[0] + " + " +
                           (self.prod2.p['filename'].split('/')[-1]).split('.')[0] + "\n\n")

        self.p['title'].append('Reactional properties in kcal/mol \n')

    def Write(self):
        for j in range(len(self.T)):
            self.p['rates'].extend(["%6s" % str("%.1f" % self.T[j]), "          ", str("%1.4E" % (1000.0 / self.T[j])),
                   "        ",
                   "%24s" % str("%1.17E" % self.p['ktst'][j]),
                   "     ", "%24s" % str("%1.17E" % self.p['kdtst'][j]), "     ", "%24s" % str("%1.17E" % self.p['kbell35'][j]),
                   "     ",
                   "%24s" % str("%1.17E" % self.p['kbell58'][j]),
                   "     ", "%24s" % str("%1.17E" % self.p['kbell2T'][j]), "     ", "%24s" % str("%1.17E" % self.p['kST'][j]),
                   " \n"])

            self.p['Lnrates'].extend(
                ["%6s" % str("%.1f" % self.T[j]), "          ", str("%1.4E" % (1000.0 / self.T[j])),
                 "        ",
                 "%24s" % str("%1.17E" % self.p['Lnktst'][j]),
                 "     ", "%24s" % str("%1.17E" % self.p['Lnkdtst'][j]), "     ", "%24s" % str("%1.17E" % self.p['Lnkbell35'][j]),
                 "     ",
                 "%24s" % str("%1.17E" % self.p['Lnkbell58'][j]), "     ", "%24s" % str("%1.17E" % self.p['Lnkbell2T'][j]),
                 "     ",
                 "%24s" % str("%1.17E" % self.p['LnkST'][j]), " \n"])
            self.p['Logrates'].extend(
                ["%6s" % str("%.1f" % self.T[j]), "          ", str("%1.4E" % (1000.0 / self.T[j])),
                 "        ",
                 "%24s" % str("%1.17E" % self.p['Logktst'][j]),
                 "     ", "%24s" % str("%1.17E" % self.p['Logkdtst'][j]), "     ", "%24s" % str("%1.17E" % self.p['Logkbell35'][j]),
                 "     ",
                 "%24s" % str("%1.17E" % self.p['Logkbell58'][j]), "     ", "%24s" % str("%1.17E" % self.p['Logkbell2T'][j]),
                 "     ",
                 "%24s" % str("%1.17E" % self.p['LogkST'][j]), " \n"])

        # ------------------------------------------------------------------------------------------------------------------------

        self.p['out'].append('       d                           Eo                Freq_imag        Tc            ΔE             ΔH             ΔG                REAC                PROD                TS\n')
        self.p['out'].extend([str("%.18f" % self.p['d']), "     ", str("%.18f" % (self.p['E0cal'] / 1000.0)), "     ", str("%.4f" % -self.p['freqi']), "     ", str("%.6f" % self.p['Tc']), "     ", str("%.6f" % self.p['de']), "     ", str("%.6f" % self.p['dh']), "     ", str("%.6f" % self.p['dg']), "     ", str("%.6f" % self.p['PES_REAC']), "     ", str("%.6f" % self.p['PES_PROD']), "     ", str("%.6f" % self.p['PES_TS']), " \n\n"])

        self.p['out'].append(" Temp              1000/T               kTST                        kdTST                       kbell35                      kbell58                      kbell582t                      kST \n")
        self.p['out'].extend(self.p['rates'])
        self.p['out'].append(
            "\n Temp             1000/T                Ln kTST                      Ln kdTST                   Ln kbell35                   Ln kbell58                      Ln kbell582t                   Ln kST\n")
        self.p['out'].extend(self.p['Lnrates'])
        self.p['out'].append(
            "\n Temp             1000/T              Log10 kTST                   Log10 kdTST                Log10 kbell35                Log10 kbell58                   Log10 kbell582t                Log10 kST\n")
        self.p['out'].extend(self.p['Logrates'])
        self.p['out'].append('\n\n')

class Smoluchowski_Calc:
    NA = 6.0221409E23
    h = 6.62607004E-34
    r = 8.3144720E0
    kb = 1.3806503E-23
    pi = (mt.atan2(1.0, 1.0)) * 4

    def __init__(self,reaction):
        self.reaction = reaction
        self.p = {}
        self.p['eta'] = [] ; self.p['Dif1'] = [] ; self.p['Dif2'] = [] ;  self.p['kD'] = []

        self.p['ktst_obs'] = [] ; self.p['kdtst_obs'] = [] ; self.p['kST_obs'] = []
        self.p['kbell35_obs'] = [] ; self.p['kbell58_obs'] = [] ; self.p['kbell2T_obs'] = []

        self.p['Lnktst_obs'] = [];self.p['Lnkdtst_obs'] = [];self.p['LnkST_obs'] = []
        self.p['Lnkbell35_obs'] = [];self.p['Lnkbell58_obs'] = [];self.p['Lnkbell2T_obs'] = []

        self.p['Logktst_obs'] = [];self.p['Logkdtst_obs'] = [];self.p['LogkST_obs'] = []
        self.p['Logkbell35_obs'] = [];self.p['Logkbell58_obs'] = [];self.p['Logkbell2T_obs'] = []


        self.p['LnDif1'] = [];self.p['LnDif2'] = [];self.p['LnkD'] = []
        self.p['LogDif1'] = [];self.p['LogDif2'] = [];self.p['LogkD'] = []

        self.p['epsilon'] = 0.0
        self.p['d'] = 0.0
        self.p['eta0'] = 0.0

        self.p['rates'] = [] ; self.p['Lnrates'] = [] ; self.p['Logrates'] = []
        self.p['out'] = []

    def Water(self):
# Fit data for water solvent
        self.p['epsilon'] = -0.58725 * self.r * 1000.0 #J/mol
        self.p['d'] = -0.3628
        self.p['eta0'] = np.exp(-8.21619) #Poise or g/cm.s

    def Calc(self):
        for j in range(len(self.reaction.T)):
            self.p['eta'].append(self.p['eta0'] * ((1 - self.p['d'] * self.p['epsilon'] / (self.r * self.reaction.T[j])) ** (1 / self.p['d'])))

            self.p['Dif1'].append((1.0E15 * self.kb * self.reaction.T[j]) / (6.0 * self.pi * self.reaction.reac1.p['rVol']) * (1.0/self.p['eta'][j]))  # cm²/s
            self.p['LnDif1'].append(np.log(self.p['Dif1'][j]))
            self.p['LogDif1'].append(np.log10(self.p['Dif1'][j]))

            self.p['Dif2'].append((1.0E15 * self.kb * self.reaction.T[j]) / (6.0 * self.pi * self.reaction.reac2.p['rVol']) * (1.0 / self.p['eta'][j]))
            self.p['LnDif2'].append(np.log(self.p['Dif2'][j]))
            self.p['LogDif2'].append(np.log10(self.p['Dif2'][j]))

            self.p['kD'].append(4.0 * self.pi * self.reaction.ts.p['rVol'] * self.NA * (self.p['Dif1'][j] + self.p['Dif2'][j]) * 1.0E-8) # cm³/mol*s
            self.p['LnkD'].append(np.log(self.p['kD'][j]))
            self.p['LogkD'].append(np.log10(self.p['kD'][j]))

            self.p['ktst_obs'].append(self.reaction.p['ktst'][j] * (self.p['kD'][j] / (self.reaction.p['ktst'][j] + self.p['kD'][j])))
            self.p['kdtst_obs'].append(self.reaction.p['kdtst'][j] * (self.p['kD'][j] / (self.reaction.p['kdtst'][j] + self.p['kD'][j])))
            self.p['kST_obs'].append(self.reaction.p['kST'][j] * (self.p['kD'][j] / (self.reaction.p['kST'][j] + self.p['kD'][j])))
            self.p['kbell35_obs'].append(self.reaction.p['kbell35'][j] * (self.p['kD'][j] / (self.reaction.p['kbell35'][j] + self.p['kD'][j])))
            self.p['kbell58_obs'].append(self.reaction.p['kbell58'][j] * (self.p['kD'][j] / (self.reaction.p['kbell58'][j] + self.p['kD'][j])))
            self.p['kbell2T_obs'].append(self.reaction.p['kbell2T'][j] * (self.p['kD'][j] / (self.reaction.p['kbell2T'][j] + self.p['kD'][j])))
            try:
                self.p['Lnktst_obs'].append(np.log(self.p['ktst_obs'][j]))
                self.p['Logktst_obs'].append(np.log10(self.p['ktst_obs'][j]))
            except:
                self.p['Lnktst_obs'].append(float("NaN"))
                self.p['Logktst_obs'].append(float("NaN"))
            try:
                self.p['Lnkdtst_obs'].append(np.log(self.p['kdtst_obs'][j]))
                self.p['Logkdtst_obs'].append(np.log10(self.p['kdtst_obs'][j]))
            except:
                self.p['Lnkdtst_obs'].append(float("NaN"))
                self.p['Logkdtst_obs'].append(float("NaN"))
            try:
                self.p['LnkST_obs'].append(np.log(self.p['kST_obs'][j]))
                self.p['LogkST_obs'].append(np.log10(self.p['kST_obs'][j]))
            except:
                self.p['LnkST_obs'].append(float("NaN"))
                self.p['LogkST_obs'].append(float("NaN"))
            try:
                self.p['Lnkbell35_obs'].append(np.log(self.p['kbell35_obs'][j]))
                self.p['Logkbell35_obs'].append(np.log10(self.p['kbell35_obs'][j]))
            except:
                self.p['Lnkbell35_obs'].append(float("NaN"))
                self.p['Logkbell35_obs'].append(float("NaN"))
            try:
                self.p['Lnkbell58_obs'].append(np.log(self.p['kbell58_obs'][j]))
                self.p['Logkbell58_obs'].append(np.log10(self.p['kbell58_obs'][j]))
            except:
                self.p['Lnkbell58_obs'].append(float("NaN"))
                self.p['Logkbell58_obs'].append(float("NaN"))
            try:
                self.p['Lnkbell2T_obs'].append(np.log(self.p['kbell2T_obs'][j]))
                self.p['Logkbell2T_obs'].append(np.log10(self.p['kbell2T_obs'][j]))
            except:
                self.p['Lnkbell2T_obs'].append(float("NaN"))
                self.p['Logkbell2T_obs'].append(float("NaN"))

    def Write(self):
        for j in range(len(self.reaction.T)):
            self.p['rates'].extend(
                ["", "%6s" % str("%.1f" % self.reaction.T[j]), "          ",
                 str("%1.4E" % (1000.0 / self.reaction.T[j])),
                 "        ",
                 "%24s" % str("%1.17E" % self.p['ktst_obs'][j]),
                 "     ", "%24s" % str("%1.17E" % self.p['kdtst_obs'][j]), "     ",
                 "%24s" % str("%1.17E" % self.p['kbell35_obs'][j]),
                 "     ",
                 "%24s" % str("%1.17E" % self.p['kbell58_obs'][j]),
                 "     ", "%24s" % str("%1.17E" % self.p['kbell2T_obs'][j]), "     ",
                 "%24s" % str("%1.17E" % self.p['kST_obs'][j]),
                 " \n"]
            )

            self.p['Lnrates'].extend(
                ["", "%6s" % str("%.1f" % self.reaction.T[j]), "          ",
                 str("%1.4E" % (1000.0 / self.reaction.T[j])),
                 "        ",
                 "%24s" % str("%1.17E" % self.p['Lnktst_obs'][j]),
                 "     ", "%24s" % str("%1.17E" % self.p['Lnkdtst_obs'][j]), "     ",
                 "%24s" % str("%1.17E" % self.p['Lnkbell35_obs'][j]),
                 "     ",
                 "%24s" % str("%1.17E" % self.p['Lnkbell58_obs'][j]),
                 "     ", "%24s" % str("%1.17E" % self.p['Lnkbell2T_obs'][j]), "     ",
                 "%24s" % str("%1.17E" % self.p['LnkST_obs'][j]),
                 " \n"]
                )
            self.p['Logrates'].extend(
                ["", "%6s" % str("%.1f" % self.reaction.T[j]), "          ",
                 str("%1.4E" % (1000.0 / self.reaction.T[j])),
                 "        ",
                 "%24s" % str("%1.17E" % self.p['Logktst_obs'][j]),
                 "     ", "%24s" % str("%1.17E" % self.p['Logkdtst_obs'][j]), "     ",
                 "%24s" % str("%1.17E" % self.p['Logkbell35_obs'][j]),
                 "     ",
                 "%24s" % str("%1.17E" % self.p['Logkbell58_obs'][j]),
                 "     ", "%24s" % str("%1.17E" % self.p['Logkbell2T_obs'][j]), "     ",
                 "%24s" % str("%1.17E" % self.p['LogkST_obs'][j]),
                 " \n"]
                )

        self.p['out'].append(' ### Smoluchowski Rates\n\n ')

        self.p['out'].append(
            "Temp       1000/T               kTSTobs                     kdTSTobs                    kbell35obs                   kbell58obs                   kbell582tobs                   kSTobs\n")
        self.p['out'].extend(self.p['rates'])
        self.p['out'].append(
            "\nTemp      1000/T                Ln kTSTobs                   Ln kdTSTobs                Ln kbell35obs                Ln kbell58obs                   Ln kbell582tobs                Ln kSTobs\n")
        self.p['out'].extend(self.p['Lnrates'])
        self.p['out'].append(
            "\nTemp      1000/T              Log10 kTSTobs                Log10 kdTSTobs             Log10 kbell35obs             Log10 kbell58obs                Log10 kbell582tobs             Log10 kSTobs\n")
        self.p['out'].extend(self.p['Logrates'])
        self.p['out'].append('\n\n')

class Kramer_Calc:
    NA = 6.0221409E23
    h = 6.62607004E-34
    r = 8.3144720E0
    kb = 1.3806503E-23
    pi = (mt.atan2(1.0, 1.0)) * 4
    c = 29979245800 #(cm/s)
    def __init__(self,reaction):
        self.reaction = reaction
        self.T = reaction.T
        self.p = {}
        self.p['eta'] = [] ; self.p['Fric'] = [] ; self.p['Kr'] = []

        self.p['ktst_obs'] = [] ; self.p['kdtst_obs'] = [] ; self.p['kST_obs'] = []
        self.p['kbell35_obs'] = [] ; self.p['kbell58_obs'] = [] ; self.p['kbell2T_obs'] = []

        self.p['Lnktst_obs'] = [];self.p['Lnkdtst_obs'] = [];self.p['LnkST_obs'] = []
        self.p['Lnkbell35_obs'] = [];self.p['Lnkbell58_obs'] = [];self.p['Lnkbell2T_obs'] = []

        self.p['Logktst_obs'] = [];self.p['Logkdtst_obs'] = [];self.p['LogkST_obs'] = []
        self.p['Logkbell35_obs'] = [];self.p['Logkbell58_obs'] = [];self.p['Logkbell2T_obs'] = []

        self.p['LnFric'] = [];self.p['LnKr'] = []
        self.p['LogFric'] = [] ; self.p['LogKr'] = []

        self.Acm = 1.0E-8 #cm
        self.mau = 1.660538E-24  #g
        self.RES = self.reaction.ts.p['rVol'] * self.Acm #cm
        self.M = self.reaction.ts.p['mass'] * self.mau    # Massa do TS em g

        self.p['epsilon'] = 0.0
        self.p['d'] = 0.0
        self.p['eta0'] = 0.0
        self.p['wFreqi'] = 2 * self.pi * self.c * self.reaction.p['freqi'] #s^-1 angular

        self.p['rates'] = [] ; self.p['Lnrates'] = [] ; self.p['Logrates'] = []
        self.p['out'] = []

    def Water(self):
#Fit data for water solvent
        self.p['epsilon']= -0.58725 * self.r * 1000.0
        self.p['d'] = -0.3628
        self.p['eta0'] = np.exp(-8.21619)

    def Calc(self):
        for j in range(len(self.T)):
            self.p['eta'].append(self.p['eta0'] * ((1 - self.p['d'] * self.p['epsilon'] / (self.r * self.T[j])) ** (1 / self.p['d'])))

            self.p['Fric'].append((6.0 * self.pi) * self.p['eta'][j] * self.RES / self.M)
            self.p['LnFric'].append(np.log(self.p['Fric'][j]))
            self.p['LogFric'].append(np.log10(self.p['Fric'][j]))

            self.p['Kr'].append((((((self.p['Fric'][j]**2.0)/4.0) + (self.p['wFreqi']**2))**(1.0/2.0)) - (self.p['Fric'][j]/2.0)) / self.p['wFreqi'])
            self.p['LnKr'].append(np.log(self.p['Kr'][j]))
            self.p['LogKr'].append(np.log10(self.p['Kr'][j]))

            self.p['ktst_obs'].append(self.p['Kr'][j] * self.reaction.p['ktst'][j])
            self.p['kdtst_obs'].append(self.p['Kr'][j] * self.reaction.p['kdtst'][j])
            self.p['kST_obs'].append(self.p['Kr'][j] * self.reaction.p['kST'][j])
            self.p['kbell35_obs'].append(self.p['Kr'][j] * self.reaction.p['kbell35'][j])
            self.p['kbell58_obs'].append(self.p['Kr'][j] * self.reaction.p['kbell58'][j])
            self.p['kbell2T_obs'].append(self.p['Kr'][j] * self.reaction.p['kbell2T'][j])

            try:
                self.p['Lnktst_obs'].append(np.log(self.p['ktst_obs'][j]))
                self.p['Logktst_obs'].append(np.log10(self.p['ktst_obs'][j]))
            except:
                self.p['Lnktst_obs'].append(float("NaN"))
                self.p['Logktst_obs'].append(float("NaN"))
            try:
                self.p['Lnkdtst_obs'].append(np.log(self.p['kdtst_obs'][j]))
                self.p['Logkdtst_obs'].append(np.log10(self.p['kdtst_obs'][j]))
            except:
                self.p['Lnkdtst_obs'].append(float("NaN"))
                self.p['Logkdtst_obs'].append(float("NaN"))
            try:
                self.p['LnkST_obs'].append(np.log(self.p['kST_obs'][j]))
                self.p['LogkST_obs'].append(np.log10(self.p['kST_obs'][j]))
            except:
                self.p['LnkST_obs'].append(float("NaN"))
                self.p['LogkST_obs'].append(float("NaN"))
            try:
                self.p['Lnkbell35_obs'].append(np.log(self.p['kbell35_obs'][j]))
                self.p['Logkbell35_obs'].append(np.log10(self.p['kbell35_obs'][j]))
            except:
                self.p['Lnkbell35_obs'].append(float("NaN"))
                self.p['Logkbell35_obs'].append(float("NaN"))
            try:
                self.p['Lnkbell58_obs'].append(np.log(self.p['kbell58_obs'][j]))
                self.p['Logkbell58_obs'].append(np.log10(self.p['kbell58_obs'][j]))
            except:
                self.p['Lnkbell58_obs'].append(float("NaN"))
                self.p['Logkbell58_obs'].append(float("NaN"))
            try:
                self.p['Lnkbell2T_obs'].append(np.log(self.p['kbell2T_obs'][j]))
                self.p['Logkbell2T_obs'].append(np.log10(self.p['kbell2T_obs'][j]))
            except:
                self.p['Lnkbell2T_obs'].append(float("NaN"))
                self.p['Logkbell2T_obs'].append(float("NaN"))

    def Write(self):
        for j in range(len(self.reaction.T)):
            self.p['rates'].extend(
                ["", "%6s" % str("%.1f" % self.reaction.T[j]), "          ",
                 str("%1.4E" % (1000.0 / self.reaction.T[j])),
                 "        ",
                 "%24s" % str("%1.17E" % self.p['ktst_obs'][j]),
                 "     ", "%24s" % str("%1.17E" % self.p['kdtst_obs'][j]), "     ",
                 "%24s" % str("%1.17E" % self.p['kbell35_obs'][j]),
                 "     ",
                 "%24s" % str("%1.17E" % self.p['kbell58_obs'][j]),
                 "     ", "%24s" % str("%1.17E" % self.p['kbell2T_obs'][j]), "     ",
                 "%24s" % str("%1.17E" % self.p['kST_obs'][j]),
                 " \n"]
            )

            self.p['Lnrates'].extend(
                ["", "%6s" % str("%.1f" % self.reaction.T[j]), "          ",
                 str("%1.4E" % (1000.0 / self.reaction.T[j])),
                 "        ",
                 "%24s" % str("%1.17E" % self.p['Lnktst_obs'][j]),
                 "     ", "%24s" % str("%1.17E" % self.p['Lnkdtst_obs'][j]), "     ",
                 "%24s" % str("%1.17E" % self.p['Lnkbell35_obs'][j]),
                 "     ",
                 "%24s" % str("%1.17E" % self.p['Lnkbell58_obs'][j]),
                 "     ", "%24s" % str("%1.17E" % self.p['Lnkbell2T_obs'][j]), "     ",
                 "%24s" % str("%1.17E" % self.p['LnkST_obs'][j]),
                 " \n"]
                )
            self.p['Logrates'].extend(
                ["", "%6s" % str("%.1f" % self.reaction.T[j]), "          ",
                 str("%1.4E" % (1000.0 / self.reaction.T[j])),
                 "        ",
                 "%24s" % str("%1.17E" % self.p['Logktst_obs'][j]),
                 "     ", "%24s" % str("%1.17E" % self.p['Logkdtst_obs'][j]), "     ",
                 "%24s" % str("%1.17E" % self.p['Logkbell35_obs'][j]),
                 "     ",
                 "%24s" % str("%1.17E" % self.p['Logkbell58_obs'][j]),
                 "     ", "%24s" % str("%1.17E" % self.p['Logkbell2T_obs'][j]), "     ",
                 "%24s" % str("%1.17E" % self.p['LogkST_obs'][j]),
                 " \n"]
                )

        self.p['out'].append(' ### Kramer Rates\n\n ')
        self.p['out'].append(
            "Temp       1000/T               kTSTobs                     kdTSTobs                    kbell35obs                   kbell58obs                   kbell582tobs                   kSTobs\n")
        self.p['out'].extend(self.p['rates'])
        self.p['out'].append(
            "\nTemp      1000/T                Ln kTSTobs                   Ln kdTSTobs                Ln kbell35obs                Ln kbell58obs                   Ln kbell582tobs                Ln kSTobs\n")
        self.p['out'].extend(self.p['Lnrates'])
        self.p['out'].append(
            "\nTemp      1000/T              Log10 kTSTobs                Log10 kdTSTobs             Log10 kbell35obs             Log10 kbell58obs                Log10 kbell582tobs             Log10 kSTobs\n")
        self.p['out'].extend(self.p['Logrates'])
        self.p['out'].append('\n\n')


########################################################################################################################




class Main:
    def __init__(self,root):

        padx = 30

        caminho = {}

### FRAMES

        framet = Frame(root)
        framet.pack(fill=BOTH, expand=True,padx=10)
        framec = LabelFrame(root,text='Espécies')
        framec.pack(fill=X, expand=True,padx=10,ipady=20)
        frameb = Frame(root)
        frameb.pack(padx=10)

        frametc = LabelFrame(framet,text='Métodos')
        frametc.pack(side=LEFT, fill=BOTH, expand=True)
        frametr = LabelFrame(framet,text='Bases')
        frametr.pack(side=LEFT, fill=BOTH, expand=True)

        framec1 = Frame(framec)
        framec1.pack(side=LEFT,fill=BOTH, expand=True)
        framec2 = Frame(framec)
        framec2.pack(side=LEFT, fill=BOTH, expand=True)
        framec3 = Frame(framec)
        framec3.pack(side=LEFT, fill=BOTH, expand=True)

        framec11 = Frame(framec1)
        framec11.pack(expand=True)
        framec12 = Frame(framec1)
        framec12.pack(expand=True)

        framec21 = Frame(framec2)
        framec21.pack(expand=True)

        framec31 = Frame(framec3)
        framec31.pack(expand=True)
        framec32 = Frame(framec3)
        framec32.pack(expand=True)

        frametct = Frame(frametc)
        frametct.pack()
        frametcc = Frame(frametc)
        frametcc.pack(fill=BOTH, expand=True)
        frametcb = Frame(frametc)
        frametcb.pack()

        frametrt = Frame(frametr)
        frametrt.pack()
        frametrc = Frame(frametr)
        frametrc.pack(fill=BOTH, expand=True)
        frametrb = Frame(frametr)
        frametrb.pack()


#### Métodos e Bases

        self.listbox_mtd = Listbox(frametcc,font=('arial', 10) ,justify=CENTER,selectmode=EXTENDED,relief='flat')
        self.listbox_mtd.pack(side=LEFT,fill=BOTH,expand=True)
        scrollbar_mtd = Scrollbar(frametcc)
        scrollbar_mtd.pack(side=RIGHT,fill=Y)
        self.listbox_mtd.config(yscrollcommand=scrollbar_mtd.set)
        scrollbar_mtd.config(command=self.listbox_mtd.yview)

        self.listbox_bs = Listbox(frametrc,font=('arial', 10) ,justify=CENTER,selectmode=EXTENDED,relief='flat')
        self.listbox_bs.pack(side=LEFT,fill=BOTH,expand=True)
        scrollbar_bs = Scrollbar(frametrc)
        scrollbar_bs.pack(side=RIGHT,fill=Y)
        self.listbox_bs.config(yscrollcommand=scrollbar_bs.set)
        scrollbar_bs.config(command=self.listbox_bs.yview)

        bt_mtd_rm = Button(frametcb,text='Remove',fg='red',command=self.mtd_rm)
        bt_mtd_rm.pack()
        bt_mtd_rm = Button(frametrb, text='Remove', fg='red',command=self.bs_rm)
        bt_mtd_rm.pack()
### Espécies
        self.DoFcp = {}
        self.DoFcp[0] = -1
        self.DoFcp[1] = -1
        self.DoFcp[2] = -1
        self.DoFcp[3] = -1
        self.DoFcp[4] = -1

        Label(framec11,text='Reagente 1').pack()
        self.ed_r1 = Entry(framec11)
        self.ed_r1.pack()
        mbr1 = Menubutton(framec11,text='Selecione o grau de liberdade',relief=RAISED)
        mbr1.menu = Menu(mbr1)
        mbr1["menu"] =  mbr1.menu
        mbr1.menu.add_command(label='Átomo',command=lambda: self.ChangeDoF(0,0))
        mbr1.menu.add_command(label='Linear',command=lambda: self.ChangeDoF(0,1))
        mbr1.menu.add_command(label='Não Linear',command=lambda: self.ChangeDoF(0,2))
        mbr1.menu.add_command(label='Extrair do arquivo', command=lambda: self.ChangeDoF(0,-1))
        mbr1.pack(pady=5)

        Label(framec12, text='Reagente 2').pack()
        self.ed_r2 = Entry(framec12)
        self.ed_r2.pack()
        mbr2 = Menubutton(framec12,text='Selecione o grau de liberdade',relief=RAISED)
        mbr2.menu = Menu(mbr2)
        mbr2["menu"] =  mbr2.menu
        mbr2.menu.add_command(label='Átomo',command=lambda: self.ChangeDoF(1,0))
        mbr2.menu.add_command(label='Linear',command=lambda: self.ChangeDoF(1,1))
        mbr2.menu.add_command(label='Não Linear',command=lambda: self.ChangeDoF(1,2))
        mbr2.menu.add_command(label='Extrair do arquivo', command=lambda: self.ChangeDoF(1,-1))
        mbr2.pack(pady=5)

        Label(framec21, text='Estado de\nTransição').pack()
        self.ed_ts = Entry(framec21)
        self.ed_ts.pack()
        mbts = Menubutton(framec21,text='Selecione o grau de liberdade',relief=RAISED)
        mbts.menu = Menu(mbts)
        mbts["menu"] =  mbts.menu
        mbts.menu.add_command(label='Átomo',command=lambda: self.ChangeDoF(2,0))
        mbts.menu.add_command(label='Linear',command=lambda: self.ChangeDoF(2,1))
        mbts.menu.add_command(label='Não Linear',command=lambda: self.ChangeDoF(2,2))
        mbts.menu.add_command(label='Extrair do arquivo', command=lambda: self.ChangeDoF(2,-1))
        mbts.pack(pady=5)

        Label(framec31, text='Produto 1').pack()
        self.ed_p1 = Entry(framec31)
        self.ed_p1.pack()
        mbp1 = Menubutton(framec31,text='Selecione o grau de liberdade',relief=RAISED)
        mbp1.menu = Menu(mbp1)
        mbp1["menu"] =  mbp1.menu
        mbp1.menu.add_command(label='Átomo',command=lambda: self.ChangeDoF(3,0))
        mbp1.menu.add_command(label='Linear',command=lambda: self.ChangeDoF(3,1))
        mbp1.menu.add_command(label='Não Linear',command=lambda: self.ChangeDoF(3,2))
        mbp1.menu.add_command(label='Extrair do arquivo', command=lambda: self.ChangeDoF(3,-1))
        mbp1.pack(pady=5)

        Label(framec32, text='Produto 2').pack()
        self.ed_p2 = Entry(framec32)
        self.ed_p2.pack()
        mbp2 = Menubutton(framec32,text='Selecione o grau de liberdade',relief=RAISED)
        mbp2.menu = Menu(mbp2)
        mbp2["menu"] =  mbp2.menu
        mbp2.menu.add_command(label='Átomo',command=lambda: self.ChangeDoF(4,0))
        mbp2.menu.add_command(label='Linear',command=lambda: self.ChangeDoF(4,1))
        mbp2.menu.add_command(label='Não Linear',command=lambda: self.ChangeDoF(4,2))
        mbp2.menu.add_command(label='Extrair do arquivo', command=lambda: self.ChangeDoF(4,-1))
        mbp2.pack(pady=5)

### Botão Selecionar Diretório e Botão Rodar

        self.listbox_dir = Listbox(frameb,font=('arial', 10) ,justify=CENTER,selectmode=EXTENDED,height=1,width=150,relief='flat')
        self.listbox_dir.pack(fill=X,expand=True)

        open = Button(frameb, text='Selecionar Diretório',font= ('Arial','12'), command=self.Open, fg='blue')
        open.pack(side=LEFT,expand=True,pady=20,padx=40)

        bt_rodar = Button(frameb,text='Rodar', font= ('Arial','12'),fg='Blue',command=self.Rodar)
        bt_rodar.pack(side=LEFT,expand=True,pady=20,padx=40)

        self.dir = ""
    def ChangeDoF(self,DoF,value):
        self.DoFcp[DoF] = value

    def mtd_rm(self):
        if len(self.listbox_mtd.curselection()) > 0:
            self.listbox_mtd.delete(self.listbox_mtd.curselection()[0],self.listbox_mtd.curselection()[-1])


    def bs_rm(self):
        if len(self.listbox_bs.curselection()) > 0:
            self.listbox_bs.delete(self.listbox_bs.curselection()[0], self.listbox_bs.curselection()[-1])

    def Open(self):
        self.dir = filedialog.askdirectory(title="Selecione o diretório onde estão as pastas com os métodos")
        if not os.path.isdir(self.dir):
            return
        self.listbox_dir.delete(0, END)
        self.listbox_dir.insert(END, self.dir)

        self.metodos = [x for x in os.listdir(self.dir) if os.path.isdir(self.dir+'/'+x)]
        self.bases = [x for x in os.listdir(self.dir+'/'+self.metodos[0]) if os.path.isdir(self.dir+'/'+self.metodos[0]+'/'+x)]

        self.listbox_mtd.delete(0,END)
        for i in range(len(self.metodos)):
            self.listbox_mtd.insert(END,self.metodos[i])

        self.listbox_bs.delete(0,END)
        for i in range(len(self.bases)):
            self.listbox_bs.insert(END, self.bases[i])


    def Rodar(self):
        if not os.path.isdir(self.dir):
            messagebox.showerror('Erro', 'Primeiramente selecione o diretório onde estão as pastas com os métodos')
        elif len(self.listbox_mtd.get(0,END)) == 0:
            messagebox.showerror('Erro', 'Selecione ou adicione ao quadro o nome das pastas dos seus métodos')
        elif len(self.listbox_bs.get(0, END)) == 0:
            messagebox.showerror('Erro', 'Selecione ou adicione ao quadro o nome das pastas das suas bases')
        elif len(self.ed_r1.get()) == 0:
            messagebox.showerror('Erro', 'Digite o nome do seu Reagente 1')
        elif len(self.ed_ts.get()) == 0:
            messagebox.showerror('Erro', 'Digite o nome do seu Estado de Transição')
        elif len(self.ed_p1.get()) == 0:
            messagebox.showerror('Erro', 'Digite o nome do seu Produto 1')
        else:
            if len(self.ed_r2.get()) == 0:
                typer2 = 'none'
            else:
                typer2 = 'reactant'
            if len(self.ed_p2.get()) == 0:
                typep2 = 'none'
            else:
                typep2 = 'product'

            species = [
                (self.ed_r1.get(),'reactant',self.DoFcp[0]),
                (self.ed_r2.get(), typer2,self.DoFcp[1]),
                (self.ed_ts.get(), 'ts',self.DoFcp[2]),
                (self.ed_p1.get(), 'product',self.DoFcp[3]),
                (self.ed_p2.get(), typep2,self.DoFcp[4])
            ]
            self.Calculate(self.dir,self.listbox_mtd.get(0,END),self.listbox_bs.get(0,END),species)


    def Calculate(self, dir, metodos, bases, species):


        s = {}  # dicionário das spécies
        r = {}  # dicionário das taxas
        solv = {}  # dicionário para taxas com efeito de solvente

        for mtd in metodos:  # loop pras pastas de métodos
            for bs in bases:  # loop pras pastas de bases
                try:
                    for sp, type, DoFcp in species:  # loop pros arquivos das espécies
                        ### Extraindo dados de cada specie ---------------------- ###
                        s[mtd + bs + sp] = Extract(dir + '/' + mtd + '/' + bs + '/' + sp, type, DoFcp)
                        if type != 'none':
                            s[mtd + bs + sp].ext()
                            if s[mtd + bs + sp].p['DoF'] == 0:
                                s[mtd + bs + sp].Atom()
                            s[mtd + bs + sp].calc_Q()
                        ## ------------------------------------------------------ ###

                    ### Calculando taxa para método e base específico ###
                    r[mtd + bs] = Rate_Calc(
                        s[mtd + bs + species[0][0]],
                        s[mtd + bs + species[1][0]],
                        s[mtd + bs + species[2][0]],
                        s[mtd + bs + species[3][0]],
                        s[mtd + bs + species[4][0]]
                    )
                    r[mtd + bs].Calc()
                    r[mtd + bs].Write()
                    ### ---------------------------------------------- ###

                    ### Calculando taxa com efeito de solvente ------------------------- ###
                    if s[mtd + bs + species[0][0]].p['rVol'] != 0.0 and s[mtd + bs + species[2][0]].p['rVol'] != 0.0 and s[mtd + bs + species[3][0]].p['rVol'] != 0.0:
                        if species[0][1] == 'reactant' and species[1][1] == 'reactant' and s[mtd + bs + species[1][0]].p['rVol'] != 0.0 and s[mtd + bs + species[4][0]].p['rVol'] != 0.0:
                            solv[mtd + bs + 'S'] = Smoluchowski_Calc(r[mtd + bs])
                            solv[mtd + bs + 'S'].Water()
                            solv[mtd + bs + 'S'].Calc()
                            solv[mtd + bs + 'S'].Write()

                        solv[mtd + bs + 'K'] = Kramer_Calc(r[mtd + bs])
                        solv[mtd + bs + 'K'].Water()
                        solv[mtd + bs + 'K'].Calc()
                        solv[mtd + bs + 'K'].Write()
                    ### ----------------------------------------------------------------  ###

                    ### Escrevendo arquivo de saída --------------------------------------###
                    out = open(self.dir + '/' + mtd + '/' + bs + '/' + 'Rates.dat', 'w', encoding="utf-8")

                    out.writelines(r[mtd + bs].p['title'])
                    out.writelines(r[mtd + bs].p['out'])
                    if species[0][1] == 'reactant' and species[1][1] == 'reactant':
                        out.writelines(solv[mtd + bs + 'S'].p['out'])
                    out.writelines(solv[mtd + bs + 'K'].p['out'])

                    out.close()
                    ### ----------------------------------------------------------------  ###
                except:
                    pass

if __name__ == '__main__':
    root = Tk()
    root.geometry('800x600+20+20')
    root.title('Transitivity Rate Auto')
    Main(root)
    root.mainloop()
