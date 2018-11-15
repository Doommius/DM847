from struct import calcsize
from struct import unpack
from struct import pack

def to_float(value):
    
    if value == 0: return 0
    val = (value & 32768) << 16;
    val |= (((value >> 10) & 31) + 112) << 23;
    val |= (value & 1023) << 13;
    s = pack('>I', val)
    return unpack('>f', s)[0]

class ims:
    def __init__(self, filename):
        self.filename = filename
        filetype = filename.split(".")[-1].lower()
        self.scale_retention = list()
        self.points = list()
        self.scale_drift = list()
        self.scale_correction = list()
        self.parameters = dict()
        #self.floatValLines = dict()
        self.ret = 0
        self.drift = 0
        self.extent = [0, 0, 0, 0]
        strPos = [0, 1, 2, 3, 4, 6, 7, 8, 13, 14, 15, 16, 17, 24, 25, 26, 30, 49, 51, 53, 55, 56, 60, 67, 75, 80, 88, 89, 99]
        floatPos  = [18, 19, 28, 29, 31, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 45, 46, 47, 50, 52, 54, 56, 57, 58, 61, 62, 63, 64, 65, 66, 68, 69, 70, 71, 72, 73, 81, 82, 83, 84, 85, 86, 87, 90, 91, 92, 93, 94, 100, 101, 102, 103, 104, 105, 107, 108, 109, 110, 111, 112, 114, 115, 116, 119]
        parameterPos = [0, 1, 2, 3, 4, 6, 7, 8, 13, 14, 15, 16, 17, 18, 19,  24, 25, 26, 28, 29, 30, 31, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 45, 46, 47, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 75, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 99, 100, 101, 102, 103, 104, 105, 107, 108, 109, 110, 111, 112, 114, 115, 116, 119]
        
        with open(filename, mode="rb") as f:
            if filetype == "csv":
            
                i = -1
                for line in f:
                    i += 1
                    linie = str(line.strip(), encoding="utf8")
                    num_lines = 0;
                    if len(linie) < 2: continue
                    qsp = linie.split(",")
                                    
                    if i == 130:
                        self.scale_retention = [float(x) for x in qsp[2:len(qsp) - 1]]
                    
                    elif i == 132:
                        self.scale_drift.append(float(qsp[0]))
                        self.scale_correction.append(float(qsp[1]))
                        # zweite Spalte: driftzeit
                        for ii in range (2, len(qsp)):
                            tmp = -float(qsp[ii])
                            #tmp = -tmp if (tmp < 0) else 0
                            lst = list()
                            lst.append(tmp)
                            self.points.append(lst)
                    
                    
                    elif i > 132:
                        self.scale_drift.append(float(qsp[0]))
                        self.scale_correction.append(float(qsp[1]))
                        
                        for ii in range (2, len(qsp)):
                            tmp = -float(qsp[ii])
                            #tmp = -tmp if (tmp < 0) else 0
                            self.points[ii - 2].append(tmp)

                    
                    else:
                        if i in parameterPos and len(qsp) > 2:
                            self.parameters[qsp[1]] = qsp[2]
                        '''
                        if floatPos.count(i) > 0 and len(qsp) > 2:
                            self.parameters[qsp[1]] float((qsp[2].split(";"))[0])
                            
                        if strPos.count(i) > 0 and len(qsp) > 2:
                            self.strValLines[qsp[1]] = qsp[2]
                        '''
                self.ret = len(self.scale_retention)
                self.drift = len(self.scale_drift)

            '''    
            elif filetype == "ims":
                data = f.read(3)
                # Länge und Strings einlesen
                for i in range(29):
                    data = f.read(2)
                    length = unpack("<h", data)[0]
                    data = f.read(length)
                    self.strValLines.append(unpack("<" + str(length) + "s", data)[0])
                    
                #Floats einlesen
                for i in range(64):
                    data = f.read(4)
                    self.floatValLines.append(unpack("<f", data)[0])
                
                data = f.read(2)
                self.ret = int(unpack("<h", data)[0])
                data = f.read(2)
                self.drift = int(unpack("<h", data)[0])
                
                
                data = f.read(1)
                bytesize = ord(unpack("<c", data)[0])
                

                # Inverse Mobilitätsskale einlesen
                for i in range(self.drift):
                    data = f.read(4)
                    self.scale_drift.append(unpack("<f", data)[0])

                # Driftzeitskala einlesen
                for i in range(self.drift):
                    data = f.read(4)
                    self.scale_correction.append(unpack("<f", data)[0])
                   
                # Punkte einlesen
                self.scale_retention = list()
                for i in range(self.ret):
                    # Retentionszeitskala einlesen
                    data = f.read(4)
                    self.scale_retention.append(unpack("<f", data)[0])
                    
                    self.points.append(list())
                    for j in range(self.drift):
                        data = f.read(bytesize)
                        if bytesize == 4:
                            a = unpack("<f", data)[0]
                            self.points[i].append(-a)
                        elif bytesize == 2:
                            a = unpack("<h", data)[0]
                            self.points[i].append(-1. * to_float(a))
            '''
            self.extent[0] = self.scale_drift[0]
            self.extent[1] = self.scale_drift[-1]
            self.extent[2] = self.scale_retention[0]
            self.extent[3] = self.scale_retention[-1]
    
    def save(self, save_filename):
        with open(save_filename, mode="wt") as sf:
            sf.write("#,data type," + self.parameters["data type"])
            for i in range(self.ret): sf.write(",")
            sf.write("\n")
            
            sf.write("#,version," + self.parameters["version"])
            for i in range(self.ret): sf.write(",")
            sf.write("\n")
            
            sf.write("#,template version," + self.parameters["template version"])
            for i in range(self.ret): sf.write(",")
            sf.write("\n")
            
            sf.write("#,AD-board type," + self.parameters["AD-board type"])
            for i in range(self.ret): sf.write(",")
            sf.write("\n")
            
            sf.write("#,ser.-no.," + self.parameters["ser.-no."])
            for i in range(self.ret): sf.write(",")
            sf.write("\n")
            
            
            sf.write("#\n")
            sf.write("#,date," + self.parameters["date"] + "\n")
            sf.write("#,time," + self.parameters["time"] + "\n")
            sf.write("#,file," + self.parameters["file"] + "\n")
            sf.write("#\n")
            sf.write("#\n")
            sf.write("#,SAMPLE INFORMATION,\n")
            sf.write("#\n")
            sf.write("#,sample type," + self.parameters["sample type"] + "\n")
            sf.write("#,sample ID," + self.parameters["sample ID"] + "\n")
            sf.write("#,comment," + self.parameters["comment"] + "\n")
            sf.write("#,location," + self.parameters["location"] + "\n")
            sf.write("#,location name," + self.parameters["location name"] + "\n")
            sf.write("#,height ASL / m," + self.parameters["height ASL / m"] + "\n")
            sf.write("#,total data acquisition time / s," + self.parameters["total data acquisition time / s"] + "\n")
            sf.write("#\n")
            sf.write("#\n")
            sf.write("#,IMS - INFORMATION,\n")
            sf.write("#\n")
            sf.write("#,operator," + self.parameters["operator"] + "\n")
            sf.write("#,operator name," + self.parameters["operator name"] + "\n")
            sf.write("#,IMS," + self.parameters["IMS"] + "\n")
            sf.write("#\n")
            sf.write("#,K0 RIP positive / cm^2/Vs," + self.parameters["K0 RIP positive / cm^2/Vs"] + "\n")
            sf.write("#,K0 RIP negative / cm^2/Vs," + self.parameters["K0 RIP negative / cm^2/Vs"] + "\n")
            sf.write("#,polarity," + self.parameters["polarity"] + "\n")
            sf.write("#,grid opening time / us," + self.parameters["grid opening time / us"] + "\n")
            sf.write("#\n")
            sf.write("#,pause / s," + self.parameters["pause / s"] + "\n")
            sf.write("#,tD interval (corr.) / ms from," + self.parameters["tD interval (corr.) / ms from"] + "\n")
            sf.write("#,tD interval (corr.) / ms to," + self.parameters["tD interval (corr.) / ms to"] + "\n")
            sf.write("#,1/K0 interval / Vs/cm^2 from," + self.parameters["1/K0 interval / Vs/cm^2 from"] + "\n")
            sf.write("#,1/K0 interval / Vs/cm^2 to," + self.parameters["1/K0 interval / Vs/cm^2 to"] + "\n")
            sf.write("#,no. of data points per spectra," + self.parameters["no. of data points per spectra"] + "\n")
            sf.write("#,no. of spectra," + self.parameters["no. of spectra"] + "\n")
            sf.write("#,no. averaged spectra," + self.parameters["no. averaged spectra"] + "\n")
            sf.write("#,baseline / signal units," + self.parameters["baseline / signal units"] + "\n")
            sf.write("#,baseline / V," + self.parameters["baseline / V"] + "\n")
            sf.write("#,V / signal unit," + self.parameters["V / signal unit"] + "\n")
            sf.write("#\n")
            sf.write("#,drift length / mm," + self.parameters["drift length / mm"] + "\n")
            sf.write("#,HV / kV," + self.parameters["HV / kV"] + "\n")
            sf.write("#,amplification / V/nA," + self.parameters["amplification / V/nA"] + "\n")
            sf.write("#\n")
            
            
            
            sf.write("#,drift gas," + self.parameters["drift gas"] + "\n")
            sf.write("#,drift gas flow / mL/min," + self.parameters["drift gas flow / mL/min"] + "\n")
            sf.write("#,sample gas," + self.parameters["sample gas"] + "\n")
            sf.write("#,sample gas flow / mL/min," + self.parameters["sample gas flow / mL/min"] + "\n")
            sf.write("#,carrier gas," + self.parameters["carrier gas"] + "\n")
            sf.write("#,carrier gas flow / mL/min," + self.parameters["carrier gas flow / mL/min"] + "\n")
            sf.write("#,pre-separation type," + self.parameters["pre-separation type"] + "\n")
            sf.write("#,pre-separation T / deg C," + self.parameters["pre-separation T / deg C"] + "\n")
            sf.write("#,sample loop T / deg C," + self.parameters["sample loop T / deg C"] + "\n")
            sf.write("#,sample loop volume / mL," + self.parameters["sample loop volume / mL"] + "\n")
            sf.write("#\n")
            
            sf.write("#,ambient T source," + self.parameters["ambient T source"] + "\n")
            sf.write("#,ambient T / deg C," + self.parameters["ambient T / deg C"] + "\n")
            sf.write("#,ambient T x^2," + self.parameters["ambient T x^2"] + "\n")
            sf.write("#,ambient T x^1," + self.parameters["ambient T x^1"] + "\n")
            sf.write("#,ambient T x^0," + self.parameters["ambient T x^0"] + "\n")
            sf.write("#,ambient T x^-1," + self.parameters["ambient T x^-1"] + "\n")
            sf.write("#,ambient T x^-2," + self.parameters["ambient T x^-2"] + "\n")
            sf.write("#,ambient p source," + self.parameters["ambient p source"] + "\n")
            sf.write("#,ambient p / hPa," + self.parameters["ambient p / hPa"] + "\n")
            sf.write("#,ambient p x^2," + self.parameters["ambient p x^2"] + "\n")
            sf.write("#,ambient p x^1," + self.parameters["ambient p x^1"] + "\n")
            sf.write("#,ambient p x^0," + self.parameters["ambient p x^0"] + "\n")
            sf.write("#,ambient p x^-1," + self.parameters["ambient p x^-1"] + "\n")
            sf.write("#,ambient p x^-2," + self.parameters["ambient p x^-2"] + "\n")
            sf.write("#\n")
            sf.write("#,6-way valve," + self.parameters["6-way valve"] + "\n")
            sf.write("#\n")
            sf.write("#\n")
            sf.write("#,EXTERNAL SAMPLING CONTROL,\n")
            sf.write("#\n")
            sf.write("#,control status," + self.parameters["control status"] + "\n")
            sf.write("#,control zero / signal units," + self.parameters["control zero / signal units"] + "\n")
            sf.write("#,control zero / V," + self.parameters["control zero / V"] + "\n")
            sf.write("#,control threshold / signal units," + self.parameters["control threshold / signal units"] + "\n")
            sf.write("#,control threshold / V," + self.parameters["control threshold / V"] + "\n")
            sf.write("#,control threshold2 / signal units," + self.parameters["control threshold2 / signal units"] + "\n")
            sf.write("#,control threshold2 / V," + self.parameters["control threshold2 / V"] + "\n")
            sf.write("#,control sampling time / s," + self.parameters["control sampling time / s"] + "\n")
            sf.write("#,control variable," + self.parameters["control variable"] + "\n")
            sf.write("#,control dimension," + self.parameters["control dimension"] + "\n")
            sf.write("#,control x^2," + self.parameters["control x^2"] + "\n")
            sf.write("#,control x^1," + self.parameters["control x^1"] + "\n")
            sf.write("#,control x^0," + self.parameters["control x^0"] + "\n")
            sf.write("#,control x^-1," + self.parameters["control x^-1"] + "\n")
            sf.write("#,control x^-2," + self.parameters["control x^-2"] + "\n")
            sf.write("#\n")
            sf.write("#\n")
            sf.write("#,STATISTICS,\n")
            sf.write("#\n")
            
            
            sf.write("#,RIP detection," + self.parameters["RIP detection"] + "\n")
            sf.write("#,tD  (RIP corr.) / ms," + self.parameters["tD  (RIP corr.) / ms"] + "\n")
            sf.write("#,1/K0 (RIP) / Vs/cm^2," + self.parameters["1/K0 (RIP) / Vs/cm^2"] + "\n")
            sf.write("#,K0 (RIP) / cm^2/Vs," + self.parameters["K0 (RIP) / cm^2/Vs"] + "\n")
            sf.write("#,SNR (RIP)," + self.parameters["SNR (RIP)"] + "\n")
            sf.write("#,WHM (RIP) / Vs/cm^2," + self.parameters["WHM (RIP) / Vs/cm^2"] + "\n")
            sf.write("#,res. power (RIP)," + self.parameters["res. power (RIP)"] + "\n")
            sf.write("#\n")
            sf.write("#,tD  (preRIP corr.) / ms," + self.parameters["tD  (preRIP corr.) / ms"] + "\n")
            sf.write("#,1/K0 (preRIP) / Vs/cm^2," + self.parameters["1/K0 (preRIP) / Vs/cm^2"] + "\n")
            sf.write("#,K0 (preRIP) / cm^2/Vs," + self.parameters["K0 (preRIP) / cm^2/Vs"] + "\n")
            sf.write("#,SNR (preRIP)," + self.parameters["SNR (preRIP)"] + "\n")
            sf.write("#,WHM (preRIP) / Vs/cm^2," + self.parameters["WHM (preRIP) / Vs/cm^2"] + "\n")
            sf.write("#,res. power (preRIP)," + self.parameters["res. power (preRIP)"] + "\n")
            sf.write("#\n")
            sf.write("#,signal RIP / V," + self.parameters["signal RIP / V"] + "\n")
            sf.write("#,signal preRIP / V," + self.parameters["signal preRIP / V"] + "\n")
            sf.write("#,RIP / preRIP," + self.parameters["RIP / preRIP"] + "\n")
            sf.write("#\n")
            sf.write("#\n")
            sf.write("#,Fims / cm^2/kV," + self.parameters["Fims / cm^2/kV"] + "\n")
            sf.write("#\n")
            sf.write("#\n")
            sf.write("#\n")
            sf.write("#\n")
            sf.write("#\n")
            sf.write("#\n")
            sf.write("#\n")
            sf.write("#\n")
            sf.write("#\n")
            sf.write("#\n")
            sf.write("\\   , tR")
            
            for i in range(self.ret): sf.write(", " + str(self.scale_retention[i]))
            sf.write("\n")
            sf.write("1/K0, tDcorr.\\SNr")
            for i in range(self.ret): sf.write(", " + str(i))
            sf.write("\n")
            
            
            for i in range(self.drift):
                sf.write(str(self.scale_drift[i]) + ", " + str(self.scale_correction[i]))
                
                for j in range(self.ret):
                    sf.write(", " + str(-self.points[j][i]))
                sf.write("\n")
            sf.write("\n")
            
class peak:
    def __init__(self):
        self.measurement_name = ""
        self.peak_name = ""
        self.r = 0
        self.t = 0
        self.index_r = 0
        self.index_t = 0
        self.signal = 0
        self.volume = 0
        self.peak_parameters = dict()
        
class peak_list:
    def __init__(self, filename=""):
        self.file_name = filename
        self.ims_peak_list = list()
        self.parameter_names = list()
        self.measurement_source = 0
        content = 0
        if filename != "":
            with open(filename, mode="rt") as f:
                content = f.read().split("\n")
            
            # check for additional parameters and store
            column_names = content[0].strip().split('\t')
            
            if len(column_names) >= 7:
                if len(column_names) > 7:
                    for i in range(7, len(column_names)):
                        self.parameter_names.append(column_names[i])
                    
                
                # reading all lines
                for i in range(1, len(content)):
                    if len(content[i]) == 0: continue
                    values = content[i].split("\t")
                    
                    if len(values) >= 7:
                        p = peak()
                        p.measurement_name = values[0]
                        p.peak_name = values[1]
                        p.t = float(values[2])
                        p.r = float(values[3])
                        p.signal = float(values[4])
                        p.volume = float(values[4])
                        p.index_t = float(values[5])
                        p.index_r = float(values[6])
                
                        # additional entries will be stored as peak parameters
                        if len(column_names) > 7:
                            for i in range(7, len(column_names)):
                                p.peak_parameters[column_names[i]] = float(values[i])
                                if column_names[i] == "volume": p.volume = float(values[i])
                        self.ims_peak_list.append(p)
                        
                    else:
                        print("Error: line", i, " in file", filename, " has incompatible format")
                        
            else:
                print("Error: file", filename, " has incompatible format")
    
    def save(self, save_filename):
        # opening output stream
        with open(save_filename, mode="wt") as f:
            # creating first row with column names
            f.write("measurement_name\tpeak_name\tt\tr\tsignal\tindex_t\tindex_r")
            if len(self.parameter_names) > 0:
                for p_name in self.parameter_names:
                    f.write("\t" + p_name)
            f.write("\n")
            
            # writing all peak values out
            for peak in self.ims_peak_list:
                f.write(peak.measurement_name + "\t")
                f.write(peak.peak_name + "\t")
                f.write(str(peak.t) + "\t")
                f.write(str(peak.r) + "\t")
                f.write(str(peak.signal) + "\t")
                f.write(str(peak.index_t) + "\t")
                f.write(str(peak.index_r))
                
                # writing if existing all additional peak parameters out
                if len(self.parameter_names) > 0:
                    for p_name in self.parameter_names:
                        f.write("\t" + str(peak.peak_parameters[p_name]))
                f.write("\n")