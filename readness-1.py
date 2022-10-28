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
import MySQLdb

class sceewv(Application):
	def __init__(self):
		Application.__init__(self, len(sys.argv), sys.argv)
		self.setDatabaseEnabled(True, True)
		self.setLoadInventoryEnabled(True)
		# self.setLoadConfigModuleEnabled(True)
		# self.setConfigModuleName("")

	def run(self):
		###1.1.1
		net = "SV"
		qcParams = ["delay","latency"]
		qcMargin = 86400
		staAcel=None
		natSts, wfIDs = self.getStas(net)
		# for wfID in wfIDs:
			# cond = 0
			# for metric in qcParams:
				# qc = self.qcquery(wfID, metric, qcMargin)
				# qcMed = np.median(qc)
				# if math.isnan(qcMed) or abs(qcMed) > 2:
					# cond = 1
			# if cond == 0:
				# if staAcel is not None:
					# staAcel["station"].append(wfID.stationCode())
					# staAcel["metric"].append(metric)
					# staAcel["value"].append(qcMed)
				# else:
					# staAcel={"station":[wfID.stationCode()],"metric":[metric], "value":[qcMed]}
		# print("National stations %s" % len(natSts))
		# print("Vertical channels of national stations with accelerometer %s" % len(wfIDs))
		# print("Number of national station with accelerometer and |delay|<2, |latency|<2: %s" % len(staAcel["station"]))
		# print(staAcel)
		# ###1.1.2
		# events, densEve = self.staDens(natSts)
		# print("%i Eventos en el catalogo " % (len(events)))
		# print("%i Eventos con 4Pth tt < 1Sth tt " % (len(densEve)))
		###1.2.1
		# staNoise = self.readJSON(["psd_value","psd_means_ratio"])
		# staBlack = self.mixSta(staNoise['station'], natSts)
		# print("National stations %s" % len(natSts))
		# print("National stations with noise alerts %s" % len(staBlack))
		# ###3.2.1
		# staTele = self.readJSON("latency")
		# staBlack = self.mixSta(staTele['station'], natSts)
		# print("National stations %s" % len(natSts))
		# print("National stations with latency alerts %s" % len(staBlack))
		# ###3.2.3
		# staQC = self.readJSON(["latency","delay","timing","gap","offset","overlap","availability","spike","rms"])
		# staBlack = self.mixSta(staQC['station'], natSts)
		# print("National stations %s" % len(natSts))
		# print("National stations with any qc alert %s" % len(staBlack))
		###4.1.1
		deltaDays=1
		today = date.today()
		endTime = today
		startTime = endTime - timedelta(days=deltaDays)
		# for sta in natSts:
			# # print("National station %s"%sta.code())
			# cond = self.ampQuery(sta.code(),startTime)
			# if cond == 1:
				# # print("Station with manual amplitude %s"%sta.code())
				# continue
			# else:
				# print("Station without manual amplitude %s"%sta.code())
		# ###4.1.2
		# # self.mvsQuery(startTime)
		# for sta in natSts:
			# # print("National station %s"%sta.code())
			# cond = self.mvsQuery(startTime,sta.code())
			# if cond == 1:
				# # print("Station with MVS amplitude %s"%sta.code())
				# continue
			# else:
				# print("Station without MVS amplitude %s"%sta.code())
		#4.2.1
		# endTime = core.Time.GMT()
		# startTime = endTime - core.TimeSpan(deltaDays*24*60*60)
		evIDs = self.IDQuery(startTime,endTime)
		self.oriQuery(evIDs)
		#4.2.3
		
		#4.2.4
		
		#4.2.5
		
		#5.1
		
		#5.2

		return True
		

	def IDQuery(self, startTime, endTime):
		evIDs = []
		eewDBuri = "mysql://sysop:sys0pm4rn@127.0.0.1:3345/seiscomp3"
		db = io.DatabaseInterface.Open(eewDBuri)
		dba = datamodel.DatabaseArchive(db)
		# build SQL query
		q = "SELECT PEvent.publicID "
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
		for i in range(IDIt.fieldCount()):
			evIDs.append(IDIt.field(i))
		try:
			while IDIt.next():
				for i in range(IDIt.fieldCount()):
					evIDs.append(IDIt.field(i))
		except:
			pass
		return evIDs

	def oriQuery(self, evIDs):
		originsID=[]
		eewDBuri = "mysql://sysop:sys0pm4rn@127.0.0.1:3345/seiscomp3"
		db = io.DatabaseInterface.Open(eewDBuri)
		query = datamodel.DatabaseQuery(db)
		for eventID in evIDs:
			cond = 0
			print("EVENT %s"%eventID)
			if cond = 0:
				for obj in query.getOriginsDescending(eventID):
					origin = datamodel.Origin.Cast(obj)
					originsID.append(origin.publicID())
				for oriID in originsID:
					origin = query.loadObject(datamodel.Origin.TypeInfo(), oriID)
					origin = datamodel.Origin.Cast(origin)
					if isinstance(origin, datamodel.Origin):
						for mag in range(origin.magnitudeCount()):
							mag = origin.magnitude(mag)
							for num_com in range(mag.commentCount()):
								comment = mag.comment(num_com)
								if comment.id() == EEW:
									print("EVENT with EEW message %s "%eventID)
									cond = 1

	def mvsQuery(self,startTime,staName):
		cond = 0
		eewDBuri = "mysql://sysop:sys0pm4rn@127.0.0.1:3345/seiscomp3"
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
		db = io.DatabaseInterface.Open(self.databaseURI())
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
		for item in stasOne:
			sta1 = item.split(".")[1]
			for station in stasSec:
				sta2 = station.code()
				if sta1 == sta2:
					if sta2 not in staBlack:
						staBlack.append(sta2)
		return staBlack

	def readJSON(self, metrics):
		jdays = np.arange(150,200,1,dtype=int)
		staNoise = None
		json_path = "/home/cam/EEW/sceewv/json_db/sceewv/scqcalert"
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

	def staDens(self, natSts):
		fdsnwsClient = "https://earthquake.usgs.gov/"
		deltaDays = 180
		minlat = 13.2
		maxlat = 14.4
		minlon = -90.2
		maxlon = -87.7
		minmag = 2.0
		maxmag = 10.0
		densEve = []
		today = date.today()
		endTime = today
		startTime = endTime - timedelta(days=deltaDays)

		client = Client(fdsnwsClient)
		cat = client.get_events(starttime=startTime,endtime=endTime,
					minmagnitude=minmag,maxmagnitude=maxmag,
					minlatitude=minlat,maxlatitude=maxlat,
					minlongitude=minlon,maxlongitude=maxlon)
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
		query = datamodel.DatabaseQuery(self.database())
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

	# def run(self):
		# self.checkMod()
		# stas = []
		# stas = [item for item in input('Enter names (net.code net2.code2) of stations installed during the last six months: ').split()]
		# print(stas)
		# now = core.Time.GMT()
		# inv = Inventory.Instance()
		# for sta in stas:
			# st = inv.getStation(sta.split(".")[0],sta.split(".")[1],now)
			# if not st:
				# print("%s is not in the inventory, please check"%sta)
			# eewSts = self.readBind(sta)
		# print("%i Stations with EEW picker binding"%len(eewSts))
		# noDat = self.staAval(eewSts)
		# print("%i Stations with EEW binding but not streaming"%len(noDat))
		# eewInvSts = []
		# means = []
		# for sta in eewSts:
			# st = inv.getStation(sta.split(".")[0],sta.split(".")[1],now)
			# eewInvSts.append(st)
			# mean = self.ppsd(sta)
			# if mean:
				# means = means+mean
		# self.staDens(eewInvSts)
		# means = np.mean(means, axis=0)
		# print(means)
		# return True

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

	# def staAval(self, eewSts):
		# noDat = []
		# client = seedClient('192.168.2.245',port=18000,timeout=1,debug=False)
		# t = UTCDateTime() - 100
		# #client = create_client('192.168.2.245:18000', on_data=self.onData)
		# #client.get_info('ID')
		# for sta in eewSts:
			# net = sta.split(".")[0]
			# code = sta.split(".")[1]
			# #print("station %s"%sta)
			# try:
				# st = client.get_waveforms(net, code, "??", "???", t, t + 5)
				# stNum = int(str(st).split()[0])
				# if stNum == 0:
					# print("No data for %s"%sta) 
					# noDat.append(sta)
			# except Exception as e:
				# print("Seedling client to %s.%s failed"%(net,code))
				# print(e)
				# continue
			# #client.select_stream(sta.split(".")[0], sta.split(".")[1], '???') 
			# #client.connect()
			# #client.run()
			# #client.close()
			# #client.select_stream("SV", "CNCH", '???')
		# return noDat
	# def onData(self, trace):
		# print('Received new data:')
		# print(trace)
		# print()
	# def staDens(self, invSts):
		# fdsnwsClient = "https://earthquake.usgs.gov/"
		# deltaDays = 180
		# minlat = 13.2
		# maxlat = 14.4
		# minlon = -90.2
		# maxlon = -87.7
		# minmag = 2.0
		# maxmag = 10.0
		# fourDist = []
		# today = date.today()
		# endTime = today
		# startTime = endTime - timedelta(days=deltaDays)

		# client = Client(fdsnwsClient)
		# cat = client.get_events(starttime=startTime,endtime=endTime,
					# minmagnitude=minmag,maxmagnitude=maxmag,
					# minlatitude=minlat,maxlatitude=maxlat,
					# minlongitude=minlon,maxlongitude=maxlon)
		# for event in cat:
			# ori = event.preferred_origin()
			# oriLat = ori.latitude
			# oriLon = ori.longitude
			# dists = False
			# for sta in invSts:
				# if sta:
					# staLat = sta.latitude()
					# staLon = sta.longitude()
					# dist = locations2degrees(oriLat, oriLon, staLat, staLon)*111.1
					# if not dists:
						# dists = {"sta":[sta.code()],"dist":[dist]}
					# else:
						# dists["sta"].append(sta.code())
						# dists["dist"].append(dist)
			# dists = pd.DataFrame(dists)
			# dists.sort_values(by='dist')
			# fourDist.append(dists.loc[3,"dist"])
		# print("4th station distances for each even%s"%fourDist)
		# print("85 percentile of 4th distances %s"%round(np.percentile(fourDist,85),2))

	# def get_miniseed(self,stream):
		# fdsnws = "http://localhost:8080" 
		# clientWf = Client(fdsnws)

		# todayUnix = int(time.time()) 
		# yesterdayUnix = todayUnix - 86400
		# today = datetime.datetime.fromtimestamp(todayUnix).strftime('%Y-%m-%d')
		# t = UTCDateTime(today)

		# yesterday = datetime.datetime.fromtimestamp(yesterdayUnix).strftime('%Y-%m-%d')
		# ty= UTCDateTime(yesterday)
		# try:
			# net = stream.split(".")[0]
			# staCode = stream.split(".")[1]
			# st = clientWf.get_waveforms(net, staCode, ",00", "HN?,EN?,SN?,HH?,EH?,SH?", ty, t)
			# inv = clientWf.get_stations(network=net, station=staCode,location=",00", channel = "HN?,EN?,SN?,HH?,EH?,SH?", level="response" )
			# print("Obtained waveform for:")
			# print(stream)
			# return st, inv
		# except Exception as e:
			# print("no waveform for:")
			# print(stream)
			# print(e)
			# return None, None
	# def ppsd(self, stream):
		# numPools = 5
		# st, inv = self.get_miniseed(stream)
		# if st == None:
			# return
		# means = []
		# for tr in st:
			# ppsd = PPSD(tr.stats, metadata=inv)
			# ppsd.add(st)
			# means.append(ppsd.get_mean()[1])
		# return means

app = sceewv()
sys.exit(app())
