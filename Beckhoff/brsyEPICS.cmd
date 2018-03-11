# This file requires a blank line at the end.

dbLoadDatabase("./tCat.dbd",0,0)
tCat_registerRecordDeviceDriver(pdbbase)
callbackSetQueueSize(10000)

tcSetScanRate(10, 5)
tcGenerateList ("C:\slowcontrols\TwinCAT3\BRS\H1BRSEY\brsy.req", "-rv -lb")
tcGenerateList ("C:\slowcontrols\TwinCAT3\BRS\H1BRSEY\brsy.ini", "-rv -l -ns")
tcLoadRecords ("C:\slowcontrols\TwinCAT3\BRS\H1BRSEY\BRS2\BRS2_Logic\BRS2_Logic.tpy", "-rv")

iocInit()
