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
		net = "SV"
		deltaDays=90
		today = date.today()
		endTime = today
		startTime = endTime - timedelta(days=deltaDays)
		natSts, wfIDs = self.getStas(net)
		###1.1.1
		if test == "1.1.1":
#			qcParams = ["delay","latency"]
			qcParams = ["latency"]
			qcMargin = 86400
			staAcel=None
			for wfID in wfIDs:
				cond = 0
				for metric in qcParams:
					qc = self.qcquery(wfID, metric, qcMargin)
					qcMed = np.median(qc)
					if math.isnan(qcMed) or abs(qcMed) > 2:
						cond = 1
				if cond == 0:
					if staAcel is not None:
						staAcel["station"].append(wfID.stationCode())
						staAcel["metric"].append(metric)
						staAcel["value"].append(qcMed)
					else:
						staAcel={"station":[wfID.stationCode()],"metric":[metric], "value":[qcMed]}
			print("National stations %s" % len(natSts))
			print("Vertical channels of national stations with accelerometer %s" % len(wfIDs))
#			print("Number of national station with accelerometer and |delay|<2, |latency|<2: %s" % len(staAcel["station"]))
			print("Number of national station with accelerometer and |latency|<2: %s" % len(staAcel["station"]))
			print(staAcel)
		###1.1.2
		elif test == "1.1.2":
			events, densEve = self.staDens(natSts)
			print("%i Eventos en el catalogo " % (len(events)))
			print("%i Eventos con 4Pth tt < 1Sth tt " % (len(densEve)))
		###1.2.1
		elif test == "1.2.1":
			staNoise = self.readJSON(["psd_value","psd_means_ratio"])
			staBlack = self.mixSta(staNoise, natSts)
			print("National stations %s" % len(natSts))
			print("National stations with noise alerts %s" % len(staBlack))
		###3.2.1
		elif test == "3.2.1":
			staTele = self.readJSON("latency")
			staBlack = self.mixSta(staTele, natSts)
			print("National stations %s" % len(natSts))
			print("National stations with latency alerts %s" % len(staBlack))
		###3.2.3
		elif test == "3.2.3":
			staQC = self.readJSON(["latency","delay","timing","gap","offset","overlap","availability","spike","rms"])
			staBlack = self.mixSta(staQC, natSts)
			print("National stations %s" % len(natSts))
			print("National stations with any qc alert %s" % len(staBlack))
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
		###5.2
		elif test == "5.2":
			self.invSlink(natSts)
		return True

	def invSlink(self, natSts):
		staSLK,staNO=[],[]
		blackLoc = ["99"]
		whiteChan = ["HN", "EN", "SN", "HH", "EH", "SH"]
		# client = seedClient('192.168.2.245',port=18000,timeout=1,debug=False)
		client = seedClient('192.168.2.245',port=18000,timeout=1,debug=False)
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
								#print("Channel not in seedlink client %s"%name)
								#print("Error: %s"%e)
								cond = 0
			if cond == 1:
				staSLK.append(sta)
			else:
				staNO.append(sta)
		print("National stations %s" % len(natSts))
		print("Stations matching with slink %s" % len(staSLK))
		print("Stations not matching with slink %s" % len(staNO))
		return

	def compEve(self,extCat,evs):
		tp,fp=[],[]
		deltaMin = 5
		deltaDist = 2
		minute_delta = timedelta(minutes=deltaMin)
		degree_delta = deltaDist
		for i in range(len(evs["time"])):
			oriTime = UTCDateTime(evs["time"][i])
			oriLat = float(evs["lat"][i])
			oriLon = float(evs["lon"][i])
			cond = 0
			for event in extCat.events:
				for origin in event.origins:
					extTime = origin.time
					deltaTime = timedelta(seconds=abs(extTime-oriTime))
					deltaLat = abs(origin.latitude - oriLat)
					deltaLon = abs(origin.longitude - oriLon)
					if deltaTime < minute_delta and deltaLat < degree_delta and deltaLon < degree_delta:
						tp.append(evs["ID"][i])
						#print("Event %s is a True Positive"%evs["ID"][i])
						cond = 1
			if cond == 0:
				fp.append(evs["ID"][i])
				#print("Event %s is not a True Positive"%evs["ID"][i])
		print("Number of events over thresholds %s"%len(evs["time"]))
		print("Number of true positives %s"%len(tp))
		print("Number of false positives %s"%len(fp))
		return

	def IDQuery(self, startTime, endTime):
		evIDs = []
		eewDBuri = "mysql://sysop:sys0pm4rn@localhost/seiscomp3"
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
		q += "AND ROUND(Origin.latitude_value,2) BETWEEN 10.0 AND 14.53 "
		q += "AND ROUND(Origin.longitude_value,2) BETWEEN -91.4 AND -86.3 "
		q += "AND Origin.depth_value BETWEEN 0 AND 80 "
		q += "AND ROUND(Magnitude.magnitude_value,1) BETWEEN 2.0 AND 10.0 "
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
		eewDBuri = "mysql://sysop:sys0pm4rn@localhost/seiscomp3"
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
									eveEEW.append(eventID)
									cond = 1
			if cond == 0:
				eveNOeew.append(eventID)
		print("Number of events over thresholds %s"%len(evIDs))
		print("Number of events with EEW message %s"%len(eveEEW))
		print("Number of events without EEW mssage %s"%len(eveNOeew))
	def mvsQuery(self,startTime,staName):
		cond = 0
		eewDBuri = "mysql://sysop:sys0pm4rn@localhost/seiscomp3"
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
		dbURIproc = "mysql://sysop:sysop@192.168.2.245/seiscomp"
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

	def mixSta(self, stasOne, stasSec):
		staBlack = []
		numAlerts = 20
		for i in range(len(stasOne["station"])):
			if stasOne["numberAlerts"][i] > numAlerts:
				sta1 = stasOne["station"][i].split(".")[1]
				for station in stasSec:
					sta2 = station.code()
					if sta1 == sta2:
						if sta2 not in staBlack:
							staBlack.append(sta2)
		return staBlack

	def readJSON(self, metrics):
		jdays = np.arange(215,305,1,dtype=int)
		staNoise = None
		json_path = "/opt/seiscomp/share/scqcalert"
		for root, dirs, files in os.walk(json_path):
			for jday in jdays:
				for item in fnmatch.filter(files, "*."+str(jday)):
					fil = os.path.join(root, item)
					with open(os.path.join(os.getcwd(), fil), 'r') as f:
						d = json.load(f)
						for alert in d["alerts"]:
							if alert["metric"] in metrics:
								if staNoise is not None:
									staNoise["station"].append(d["mseedid"])
									staNoise["metric"].append(alert["metric"])
									staNoise["numberAlerts"].append(len(alert["values"]))
								else:
									staNoise = {"station":[d["mseedid"]],"metric":[alert["metric"]],"numberAlerts":[len(alert["values"])]}
		return staNoise

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
		deltaDays = 90
		minlat = 9.5
		maxlat = 15.0
		minlon = -92.0
		maxlon = -86.0
		minmag = 2.0
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
		densEve = []
		cat = self.extCat()
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
			dists.sort_values(by='dist')
			dist1S = dists.loc[0,"dist"]
			dist4P = dists.loc[3,"dist"]
			pFourTime,sFstTime = self.tavTimes(oriDep,dist4P,dist1S)
			if pFourTime < sFstTime:
				densEve.append(event)
		return cat, densEve
 
	def qcquery(self, stName, param, qcMargin):
		qc_vec = []
		startTime = core.Time.GMT() - core.TimeSpan(qcMargin)
		endTime = core.Time.GMT()
		dbURIproc = "mysql://sysop:sysop@192.168.2.245/seiscomp"
		db = io.DatabaseInterface.Open(dbURIproc)
		#dba = datamodel.DatabaseArchive(db)
		query = datamodel.DatabaseQuery(db)
		for obj in query.getWaveformQuality(stName,param,startTime,endTime):
			qc = datamodel.WaveformQuality.Cast(obj)
			qc_vec.append(qc.value())
		return qc_vec

	def getStas(self, netCode):
		now = core.Time.GMT()
		stations, stasIDs, wfIDs = [], [], []
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
			wfID = datamodel.WaveformStreamID()
			wfID.setNetworkCode(sta.split(".")[0])
			wfID.setStationCode(sta.split(".")[1])
			wfID.setLocationCode(sta.split(".")[2])
			wfID.setChannelCode(sta.split(".")[3])
			wfIDs.append(wfID)
		return stations, wfIDs

	# def checkMod(self):
		# ei = system.Environment.Instance()
		# root = ei.installDir()
		# shouldRun = True
		# modules = ["scmaster", "fdsnws", "scevent", "scimport", "scamp", "scmag",
			# "scvsnpick", "scvspick", "scvsloc", "scvsnloc", "scanloc", 
			# "screlocan", "screloc", "sceewlog", "sccaplog", "sceewmail", 
			# "scfdeewd", "sceewenv", "scvsmag", "scfd20asym", "scfd85sym", 
			# "scfdcrust", "scqcalert"]
		# env = Environment(root)
		# for module in modules:
			# md = Module(env,module)
			# enbl = env.isModuleEnabled(module)
			# run = md.isRunning()
			# # md.status(shouldRun)
			# if not enbl:
				# print("%s not enabled, please check" % module)
			# if not run:
				# print("%s not running please check" % module)

	# def readBind(self, sta):
		# eewMods =  ["scvspick", "scvsnpick"]
		# cfgMod = None
		# config = client.ConfigDB.Instance().config()
		# moduleName = self.name()
		# datamodel.Notifier.Enable()
		# cond = False
		# eewSts = []
		# for mod in eewMods:
			# for i in range(config.configModuleCount()):
				# cfgMod = config.configModule(i)
				# for j in range(cfgMod.configStationCount()):
					# cs = cfgMod.configStation(j)
					# for k in range(cs.setupCount()):
						# stp = cs.setup(k)
						# if stp.name() == mod and stp.enabled() == True:
							# name = "%s.%s"%(cs.networkCode(),cs.stationCode())
							# if name not in eewSts:
								# eewSts.append(name)
							# if name == sta:
								# cond = True
								# print("%s bind with %s"%(sta,mod))
								# paramSet = datamodel.ParameterSet.Find(stp.parameterSetID())
								# for l in range(paramSet.parameterCount()):
									# p_i= paramSet.parameter(l)
# #									print(p_i.name())
# #									print(p_i.value())
			# if not cond:
				# print("%s not bind with %s, please check"%(sta,mod))
		# return eewSts

app = sceewv()
sys.exit(app())
