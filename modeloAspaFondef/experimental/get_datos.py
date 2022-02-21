from subprocess import call
import glob

def get_datos(fname, datadir="./datos"):
	files = glob.glob(datadir + "/" + fname)
	if len(files) > 0:
		print ("Found ", fname, " in datadir: ", datadir)
	else:
		print ("Getting  ", fname, " from data repository. gdrive line:")
		callstring = "/home/susa/bin/gdrive download query \"name contains '"+fname+"'\" --path ./datos"
		print (callstring)
		call([callstring],shell=True)
