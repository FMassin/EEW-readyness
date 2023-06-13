from obspy.clients.seedlink.easyseedlink import create_client
from obspy.clients.seedlink.basic_client import Client as seedClient
from datetime import date, timedelta
from obspy import UTCDateTime
from obspy.clients.fdsn import Client
from obspy.geodetics.base import locations2degrees
from obspy.signal import PPSD
from seiscomp import client, core, datamodel, shell, config, system, io
from seiscomp.client import Inventory, Application
from seiscomp.kernel import Environment, Module
from multiprocessing import Pool
from obspy.taup import TauPyModel
import pandas as pd
import numpy as np
import sys, os, datetime, time, math, fnmatch, json

class sceewv(Application):
	def __init__(self):
		Application.__init__(self, len(sys.argv), sys.argv)
		self.setDatabaseEnabled(True, True)
		self.setLoadInventoryEnabled(True)
		# self.setLoadConfigModuleEnabled(True)
		# self.setConfigModuleName("")

	def run(self):
                test = "5.2"
                net = "GI"
                deltaDays=180
                today = date.today()
                endTime = today
                startTime = endTime - timedelta(days=deltaDays)
                natSts, wfIDs, staNoRep = self.getStas(net)
                ###1.1.1
                if test == "1.1.1":
#                        qcParams = ["delay","latency"]
                        qcParams = ["delay"]
                        qcMargin = 86400
                        staAcel,staAcelNo,staNAN=None,None,None
                        staAcelNAN,staNAN = [],[]
                        for wfID in wfIDs:
                                cond = 0
                                for metric in qcParams:
                                        qc = self.qcquery(wfID, metric, qcMargin)
                                        qcMed = np.median(qc)
                                        if math.isnan(qcMed):
                                                cond = 1
                                        if abs(qcMed) > 2:
                                                cond = 2
                                if cond == 0:
                                        if staAcel is not None:
                                                staAcel["station"].append(wfID.stationCode())
                                                staAcel["metric"].append(metric)
                                                staAcel["value"].append(qcMed)
                                        else:
                                                staAcel={"station":[wfID.stationCode()],"metric":[metric], "value":[qcMed]}
                                if cond == 1:
                                        staNAN.append(wfID.stationCode())
                                if cond == 2:
                                        if staAcelNo is not None:
                                                staAcelNo["station"].append(wfID.stationCode())
                                                staAcelNo["metric"].append(metric)
                                                staAcelNo["value"].append(qcMed)
                                        else:
                                                staAcelNo={"station":[wfID.stationCode()],"metric":[metric], "value":[qcMed]}
                        for sta in staNAN:
                                if sta not in staAcelNo["station"] and sta not in staAcel["station"] and sta not in staAcelNAN:
                                        staAcelNAN.append(sta)

                        print("National stations %s" % len(natSts))
                        print("National stations with accelerometer %s" % len(staNoRep))
                        print(staNoRep)
#                       print("Number of national station with accelerometer and |delay|<2, |latency|<2: %s" % len(staAcel["station"]))
                        print("Number of national station with accelerometer and |delay|<2: %s" % len(staAcel["station"]))
                        print(staAcel)
                        print("Number of national station with accelerometer and |delay|>2: %s" % len(staAcelNo["station"]))
                        print(staAcelNo)
                        print("Estaciones sin dataos de qc %s"% len(staAcelNAN))
                        print(staAcelNAN)
                ###1.1.2
                elif test == "1.1.2":
                        events, densEve = self.staDens(natSts)
                        print(events)
                        print("%i Eventos en el catalogo " % (len(events)))
                        print("%i Eventos con 4Pth tt < 1Sth tt " % (len(densEve)))
                ###1.2.1
                elif test == "1.2.1":
                        numJdays = 30
                        staNoise = self.readJSON(["psd_value","psd_means_ratio"],numJdays)
                        numAlerts = 5
                        psdSts = self.readSTAPSD()
                        staBlack = self.mixSta(staNoise, natSts, numAlerts, numJdays)
                        staMiss = self.misSta(psdSts, natSts)
                        print("National stations with no PSDs computation %s" % len(staMiss))
                        print(staMiss)
                        print("National stations %s" % len(natSts))
                        print("National stations with noise alerts %s" % len(staBlack))
                        print(staBlack)
                ###3.2.1
                elif test == "3.2.1":
                        numJdays = 30
                        staTele = self.readJSON("latency",numJdays)
                        numAlerts = 20
                        staBlack = self.mixSta(staTele, natSts, numAlerts, numJdays)
                        print("National stations %s" % len(natSts))
                        print("National stations with latency alerts %s" % len(staBlack))
                        print(staBlack)
                ###3.2.3
                elif test == "3.2.3":
                        numJdays = 30
                        numAlerts = 20
                        staQC = self.readJSON(["latency","delay","timing","gap","offset","overlap","availability","spike","rms"], numJdays)
                        staBlack = self.mixSta(staQC, natSts, numAlerts, numJdays)
                        print("National stations %s" % len(natSts))
                        print("National stations with any qc alert %s" % len(staBlack))
                        print(staBlack)
                ###4.1.1
                elif test == "4.1.1":
                        staMAN,staNO = [],[]
                        for sta in natSts:
                                # print("National station %s"%sta.code())
                                cond = self.ampQuery(sta.code(),startTime)
                                if cond == 1:
                                        staMAN.append(sta.code())
                                        #print("Station with manual amplitude %s"%sta.code())
                                        #continue
                                else:
                                        staNO.append(sta.code())
                                        #print("Station without manual amplitude %s"%sta.code())
                        print("National stations %s" % len(natSts))
                        print("National stations with manual amplitudes %s" % len(staMAN))
                        print("National stations without manual amplitudes %s" % len(staNO))
                        print(staMAN)
                        print(staNO)
                ###4.1.2
                elif test == "4.1.2":
                        staMVS,staNO=[],[]
                        for sta in natSts:
                                # print("National station %s"%sta.code())
                                cond = self.mvsQuery(startTime,sta.code())
                                if cond == 1:
                                        staMVS.append(sta.code())
                                        #print("Station with MVS amplitude %s"%sta.code())
                                        #continue
                                else:
                                        staNO.append(sta.code())
                                        #print("Station without MVS amplitude %s"%sta.code())
                        print("National stations %s" % len(natSts))
                        print("National stations with MVS amplitudes %s" % len(staMVS))
                        print("National stations without MVS amplitudes %s" % len(staNO))
                        print(staMVS)
                        print(staNO)
                ##4.2.1
                elif test == "4.2.1":
                        #endTime = core.Time.GMT()
                        #startTime = endTime - core.TimeSpan(deltaDays*24*60*60)
                        evs = self.IDQuery(startTime,endTime)
                        self.oriQuery(evs["ID"])
                ##4.2.2
                elif test == "4.2.2":
                        evs = self.IDQuery(startTime,endTime)
                        extCat = self.extCat()
                        self.compEve(extCat,evs)
                ###4.2.3
                elif test == "4.2.3":
                        evs = self.IDQuery(startTime,endTime)
                        extCat = self.extCat()
                        self.compEve(extCat,evs)
                ###4.2.4
                ###5.1
                elif test == "5.1":
                        evsMAN,evsAUT=[],[]
                        evs = self.IDQuery(startTime,endTime)
                        for i in range(len(evs["ID"])):
                                if evs["evalMode"][i] == "manual":
                                        evsMAN.append(evs["ID"][i])
                                        #print("Event %s with preferred manual solution"%evs["ID"][i])
                                else:
                                        evsAUT.append(evs["ID"][i])
                                        #print("Event %s without preferred manual solution"%evs["ID"][i])
                        print("Events over thresholds %s" % len(evs["ID"]))
                        print("Events with preferred manual solution %s" % len(evsMAN))
                        print("Events with preferred automatic solution %s" % len(evsAUT))
                        print(evsMAN)
                        print(evsAUT)
                ###5.2
                elif test == "5.2":
                        self.invSlink(natSts)
                return True

	def invSlink(self, natSts):
		staSLK,staNO=[],[]
		blackLoc = ["99"]
		whiteChan = ["HN", "EN", "SN", "HH", "EH", "SH"]
		# client = seedClient('192.168.2.245',port=18000,timeout=1,debug=False)
		client = seedClient('172.20.9.20',port=18000,timeout=1,debug=False)
		t = UTCDateTime() - 100
		for station in natSts:
			cond = 1
			sta = station.code()
			locs = station.sensorLocationCount()
			for l in range(locs):
				location = station.sensorLocation(l)
				loc = location.code()
				if loc not in blackLoc:
					chans = location.streamCount()
					for c in range(chans):
						channel = location.stream(c)
						chan = channel.code()
						if chan[0:2] in whiteChan:
							name=[sta,loc,chan]
							try:
								info = client.get_info(station=sta,location=loc,channel=chan,level="channel")
								if len(info) > 0:
									print("Channel in seedlink client %s"%info)
								else:
									#print("Channel not in seedlink client %s"%name)
									cond = 0
							except Exception as e:
								print("Channel not in seedlink client %s"%name)
								print("Error: %s"%e)
								cond = 0
			if cond == 1:
				staSLK.append(sta)
			else:
				staNO.append(sta)
		print("National stations %s" % len(natSts))
		print("Stations matching with slink %s" % len(staSLK))
		print("Stations not matching with slink %s" % len(staNO))
		print(staSLK)
		print(staNO)
		return

	def compEve(self,extCat,evs):
                tp,fp,missed=[],[],[]
                seconds = 300
                deltaDist = 1.5
                deltaSec = timedelta(seconds=seconds)
                degree_delta = deltaDist
                for event in extCat.events:
                        cond = 0
                        for i in range(len(evs["time"])):
                                oriTime = UTCDateTime(evs["time"][i])
                                oriLat = float(evs["lat"][i])
                                oriLon = float(evs["lon"][i])
                                for origin in event.origins:
                                        extTime = origin.time
                                        deltaTime = timedelta(seconds=abs(extTime-oriTime))
                                        deltaLat = abs(origin.latitude - oriLat)
                                        deltaLon = abs(origin.longitude - oriLon)
                                        if deltaTime < deltaSec and deltaLat < degree_delta and deltaLon < degree_delta:
                                                if evs["ID"][i] not in tp:
                                                         cond = 1
                                                         tp.append(evs["ID"][i])
                                        if cond == 1:
                                                break
                                else:
                                        continue
                                break
                                                #print("Event %s is a True Positive"%evs["ID"][i])
                        if cond == 0:
                                missed.append(event.preferred_origin().time)
        #               if cond == 0:
        #                       fp.append(evs["ID"][i])
                                #print("Event %s is not a True Positive"%evs["ID"][i])
                print("Number of events EEW catalog %s"%len(evs["time"]))
                print("Number of events external catalog %s"%len(extCat.events))
                print("Number of true positives %s"%len(tp))
        #       print("Number of false positives %s"%len(fp))
                print("Missed events %s"%len(missed))
                print(tp)
                print(missed)
                return

	def IDQuery(self, startTime, endTime):
		evIDs = []
		eewDBuri = "mysql://sysop:sysop@localhost/seiscomp"
		db = io.DatabaseInterface.Open(eewDBuri)
		dba = datamodel.DatabaseArchive(db)
		# build SQL query
		q = "SELECT PEvent.publicID, Origin.latitude_value, Origin.longitude_value, Origin.depth_value, Magnitude.magnitude_value, Origin.time_value, Origin.evaluationMode "
		q += "FROM Event, PublicObject AS PEvent, Origin, PublicObject AS POrigin, Magnitude, PublicObject AS PMagnitude "
		q += "WHERE Origin._oid = POrigin._oid "
		q += "AND Event._oid = PEvent._oid "
		q += "AND Magnitude._oid = PMagnitude._oid "
		q += "AND Magnitude._parent_oid = Origin._oid "
		q += "AND Event.preferredOriginID = POrigin.publicID "
		q += "AND Event.preferredMagnitudeID = PMagnitude.publicID "
		q += "AND Origin.time_value BETWEEN '{0}' AND '{1}' "\
			.format(startTime,endTime)
		q += "AND ROUND(Origin.latitude_value,2) BETWEEN 12.2 AND 18.3 "
		q += "AND ROUND(Origin.longitude_value,2) BETWEEN -93.0 AND -89.2 "
		q += "AND ROUND(Magnitude.magnitude_value,1) BETWEEN 4.0 AND 10.0 "
		IDIt = dba.getObjectIterator(q, datamodel.Event.TypeInfo())
		evs={"ID":[IDIt.field(0)],"lat":[IDIt.field(1)],"lon":[IDIt.field(2)],"deep":[IDIt.field(3)],"mag":[IDIt.field(4)],"time":[IDIt.field(5)],"evalMode":[IDIt.field(6)]}
		try:
			while IDIt.next():
				for i in range(IDIt.fieldCount()):
					if i == 0:
						evs["ID"].append(IDIt.field(i))
					elif i == 1:
						evs["lat"].append(IDIt.field(i))
					elif i == 2:
						evs["lon"].append(IDIt.field(i))
					elif i == 3:
						evs["deep"].append(IDIt.field(i))
					elif i == 4:
						evs["mag"].append(IDIt.field(i))
					elif i == 5:
						evs["time"].append(IDIt.field(i))
					elif i == 6:
						evs["evalMode"].append(IDIt.field(i))
		except:
			pass
		return evs

	def oriQuery(self, evIDs):
                eveEEW,eveNOeew=[],[]
                eewDBuri = "mysql://sysop:sysop@localhost/seiscomp"
                db = io.DatabaseInterface.Open(eewDBuri)
                query = datamodel.DatabaseQuery(db)
                for eventID in evIDs:
                        originsID=[]
                        cond = 0
                        #print("EVENT %s"%eventID)
                        for obj in query.getOriginsDescending(eventID):
                                origin = datamodel.Origin.Cast(obj)
                                originsID.append(origin.publicID())
                        for oriID in originsID:
                                if cond == 0:
                                        origin = query.loadObject(datamodel.Origin.TypeInfo(), oriID)
                                        origin = datamodel.Origin.Cast(origin)
                                        if isinstance(origin, datamodel.Origin):
                                                for mag in range(origin.magnitudeCount()):
                                                        mag = origin.magnitude(mag)
                                                        for num_com in range(mag.commentCount()):
                                                                comment = mag.comment(num_com)
                                                                if comment.id() == 'EEW':
                                                                        #print("EVENT with EEW message %s "%eventID)
                                                                        #print(comment.id(),comment.text())
                                                                        cond = 1
                                                                        if eventID not in eveEEW:
                                                                                eveEEW.append(eventID)
                        if cond == 0:
                                eveNOeew.append(eventID)
                print("Number of events over thresholds %s"%len(evIDs))
                print("Number of events with EEW message %s"%len(eveEEW))
                print("Number of events without EEW mssage %s"%len(eveNOeew))
                print(eveEEW)
                print(eveNOeew)
	def mvsQuery(self,startTime,staName):
		cond = 0
		eewDBuri = "mysql://sysop:sysop@localhost/seiscomp"
		db = io.DatabaseInterface.Open(eewDBuri)
		dba = datamodel.DatabaseArchive(db)
		# build SQL query
		q = "SELECT * from StationMagnitude " \
			"WHERE type = 'MVS' "
		q +="AND creationInfo_creationTime >= '{0}' " \
			.format(startTime)
		q += "AND waveformID_stationCode = '{0}' "\
			.format(staName)
		ampIt = dba.getObjectIterator(q, datamodel.Amplitude.TypeInfo())
		amp = datamodel.Amplitude.Cast(ampIt.get())
		while cond == 0:
			if amp is None:
				break
			else:
				cond = 1
				return cond

	def ampQuery(self,staName,startTime):
		cond = 0
		evalMode = "0, manual"
		dbURIproc = "mysql://sysop:sysop@172.20.11.5/seiscomp"
		db = io.DatabaseInterface.Open(dbURIproc)
		dba = datamodel.DatabaseArchive(db)
		# build SQL query
		q = "SELECT * from Amplitude " \
			"WHERE amplitude_used = '1' "
		q +="AND creationInfo_creationTime >= '{0}' " \
			.format(startTime)
		q += "AND evaluationMode = 'manual' "
		q += "AND waveformID_stationCode = '{0}' "\
			.format(staName)
		ampIt = dba.getObjectIterator(q, datamodel.Amplitude.TypeInfo())
		amp = datamodel.Amplitude.Cast(ampIt.get())
		while cond == 0:
			if amp is None:
				break
			else:
				cond = 1
				return cond

	def mixSta(self, stasOne, stasSec, numAlerts, numJdays):
                staBlack = []
                dayRatio = 0.1
                for sta in stasOne:
                        days = []
                        for i in range(len(stasOne[sta])):
                                if stasOne[sta][i]["numberAlerts"] > numAlerts:
                                        day = stasOne[sta][i]["day"]
                                        if day not in days:
                                                days.append(day)
                        if len(days) != 0 and len(days)/numJdays > dayRatio:
                                staBlack.append(sta)
#                               sta1 = stasOne["station"][i].split(".")[1]
#                               for station in stasSec:
#                                       sta2 = station.code()
#                                       if sta1 == sta2:
#                                               if sta2 not in staBlack:
#                                                       staBlack.append(sta2)
                return staBlack
	def misSta(self, stasOne, stasSec):
                misSts, sta1Vec = [], []
                for station in stasSec:
                        cond = 0
                        sta2 = station.code()
                        for sta1 in stasOne:
                                if sta1 == sta2:
                                        cond = 1
                                        if sta1 in sta1Vec:
                                                print(sta1)
                                        sta1Vec.append(sta1)
                        if cond == 0:
                                misSts.append(sta2)
                return misSts

	def readJSON(self, metrics, numJdays):
                startDay = 1
                endDay = startDay+numJdays
                jdays = np.arange(startDay,endDay,1,dtype=int)
                staNoise = {}
                staVec, staVecNoise = [], []
                json_path = "/opt/seiscomp/share/scqcalert/2023/GI"
                for root, dirs, files in os.walk(json_path):
                        for jday in jdays:
                                for item in fnmatch.filter(files, "*.2023."+str(jday)):
                                        fil = os.path.join(root, item)
                                        with open(os.path.join(os.getcwd(), fil), 'r') as f:
                                                filedata = f.read()
                                                if '\n}{\n' in filedata:
                                                        filedata = filedata.replace('\n}{\n', ',\n')
                                                        with open(os.path.join(os.getcwd(), fil), 'w') as f:
                                                                f.write(filedata)
                                        with open(os.path.join(os.getcwd(), fil), 'r') as f:
                                                d = json.load(f)
                                                sta = d["mseedid"].split(".")[1]
                                                if sta not in staVec:
                                                        staVec.append(sta)
                                                for alert in d["alerts"]:
                                                        if alert["metric"] in metrics:
                                                                if sta not in staVecNoise:
                                                                        staVecNoise.append(sta)
                                                                staDict = {"metric":alert["metric"], "numberAlerts":len(alert["values"]), "day":jday}
                                                                if sta not in staNoise:
                                                                        staNoise[sta]=[staDict]
                                                                else:
                                                                        staNoise[sta].append(staDict)
                print("Estaciones con alertas de scqclaert")
                print(staVec)
                print(len(staVec))
                print("Estaciones con alertas en metrica evaluada de scqclaert")
                print(staVecNoise)
                print(len(staVecNoise))
                return staNoise

	def readSTAPSD(self):
                staVec=[]
                jdays = np.arange(1,26,1,dtype=int)
                staNoise = None
                json_path = "/home/insivumeh/PSD/PNG/2023/GI"
                for root, dirs, files in os.walk(json_path):
                        for jday in jdays:
                                jday = str(jday).zfill(3)
                                for item in fnmatch.filter(files, "*.2023."+jday+".npz"):
                                        sta=item.split(".")[1]
                                        if sta not in staVec:
                                                staVec.append(sta)
                print("Estaciones con PSDs")
                print(staVec)
                print(len(staVec))
                return staVec

	def tavTimes(self, oriDep, dist4P, dist1S):
		MODEL = 'iasp91'
		model = TauPyModel(model=MODEL)
		pArr = model.get_travel_times(source_depth_in_km=oriDep,distance_in_degree=dist4P, phase_list=["p", "P"])
		sArr = model.get_travel_times(source_depth_in_km=oriDep,distance_in_degree=dist1S, phase_list=["s", "S"])
		pFourTime = float(str(pArr[0]).split(" ")[4])
		sFstTime = float(str(sArr[0]).split(" ")[4])
		return pFourTime,sFstTime

	def extCat(self):
		fdsnwsClient = "https://earthquake.usgs.gov/"
		deltaDays = 180
		minlat = 12.4
		maxlat = 18.1
		minlon = -92.8
		maxlon = -89.4
		minmag = 4.0
		maxmag = 10.0
		today = date.today()
		endTime = today
		startTime = endTime - timedelta(days=deltaDays)
		client = Client(fdsnwsClient)
		cat = client.get_events(starttime=startTime,endtime=endTime,
					minmagnitude=minmag,maxmagnitude=maxmag,
					minlatitude=minlat,maxlatitude=maxlat,
					minlongitude=minlon,maxlongitude=maxlon)
		return cat

	def staDens(self, natSts):
                densEve,pFourTimeVec,sFstTimeVec = [],[],[]
                cat = self.extCat()
                #cat.plot(projection="local",outfile="exCat.png")
                for event in cat:
                        ori = event.preferred_origin()
                        oriLat = ori.latitude
                        oriLon = ori.longitude
                        oriDep = ori.depth/1000
                        dists = False
                        for sta in natSts:
                                if sta:
                                        staLat = sta.latitude()
                                        staLon = sta.longitude()
                                        dist = locations2degrees(oriLat, oriLon, staLat, staLon)
                                        if not dists:
                                                dists = {"sta":[sta.code()],"dist":[dist]}
                                        else:
                                                dists["sta"].append(sta.code())
                                                dists["dist"].append(dist)
                        dists = pd.DataFrame(dists)
                        distSort = dists.sort_values(by='dist',ignore_index=True)
                        dist1S = distSort.loc[0,"dist"]
                        dist4P = distSort.loc[3,"dist"]
                        pFourTime,sFstTime = self.tavTimes(oriDep,dist4P,dist1S)
                        pFourTimeVec.append(pFourTime)
                        sFstTimeVec.append(sFstTime)
                        if pFourTime < sFstTime:
                                densEve.append(event)
                print(pFourTimeVec)
                print(sFstTimeVec)
                return cat, densEve
 
	def qcquery(self, stName, param, qcMargin):
		qc_vec = []
		startTime = core.Time.GMT() - core.TimeSpan(qcMargin)
#		startTime = core.Time.GMT() - core.TimeSpan(86400*60) - core.TimeSpan(qcMargin)
		endTime = core.Time.GMT()
#		dbURIproc = "mysql://sysop:sysop@172.20.9.20/seiscomp3"
		dbURIproc = "mysql://sysop:sysop@172.20.11.5/seiscomp"
		db = io.DatabaseInterface.Open(dbURIproc)
		#dba = datamodel.DatabaseArchive(db)
		query = datamodel.DatabaseQuery(db)
		for obj in query.getWaveformQuality(stName,param,startTime,endTime):
			qc = datamodel.WaveformQuality.Cast(obj)
			qc_vec.append(qc.value())
		return qc_vec

	def getStas(self, netCode):
                now = core.Time.GMT()
                stations, stasIDs, wfIDs, staNoRep = [], [], [], []
                blackLoc = ["99"]
                acelChan = ["HNZ", "ENZ", "SNZ"]
                inv = Inventory.Instance()
                nets = inv.inventory().networkCount()
                for net in range(nets):
                        network = inv.inventory().network(net)
                        if network.code() == netCode:
                                stas = network.stationCount()
                                for s in range(stas):
                                        station = network.station(s)
                                        sta = station.code()
                                        station = inv.getStation(netCode,sta,now)
                                        if station and station not in stations:
                                                stations.append(station)
                for station in stations:
                        sta = station.code()
                        locs = station.sensorLocationCount()
                        for l in range(locs):
                                location = station.sensorLocation(l)
                                loc = location.code()
                                if loc not in blackLoc:
                                        chans = location.streamCount()
                                        for c in range(chans):
                                                channel = location.stream(c)
                                                chan = channel.code()
                                                if chan in acelChan:
                                                        staID = "%s.%s.%s.%s" % (netCode,sta,loc,chan)
                                                        if staID not in stasIDs:
                                                                stasIDs.append(staID)
                for sta in stasIDs:
                        name = sta.split(".")[1]
                        if name not in staNoRep:
                                staNoRep.append(name)
                        wfID = datamodel.WaveformStreamID()
                        wfID.setNetworkCode(sta.split(".")[0])
                        wfID.setStationCode(sta.split(".")[1])
                        wfID.setLocationCode(sta.split(".")[2])
                        wfID.setChannelCode(sta.split(".")[3])
                        wfIDs.append(wfID)
                return stations, wfIDs, staNoRep
app = sceewv()
sys.exit(app())


