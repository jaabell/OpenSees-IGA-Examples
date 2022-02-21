
#Filtering pass-band [f1,  f2]
f1 = 2.
# f1 = 0.5
f2 = 125

f1=0.15
f2 = 140


# filtet_type = "sos" (filtrar) o "none" (no filtrar)

#["1901021_12_55_59.txt",     ( 0., 10.),     "none" ,         "view"      ]
 
# cases = [
# #     filename.txt              (t1, t2 )    filter_type        action             COMMENT
# #                                                        view|save|skip|todo|done
#     ["190713_18_07_50.txt",     (1., 16.),     "sos" ,         "done"      ],  #Blade 7:Fatigue test, 41 hours at 4.48 Hz pt.1
#     ["190715_13_20_05.txt",     (19., 29.),     "sos" ,         "save"      ],  #Blade 7:Fatigue test, 41 hours at 4.48 Hz pt.1
    
# ]


cases = [
#     filename.txt              (t1, t2 )    filter_type        action             COMMENT
#                                                        view|save|skip|todo|done
    # ["190705_18_29_09.txt",     (11.9, 30.),     "sos" ,         "save"      ],  #Blade 7:Fatigue test, 41 hours at 4.48 Hz pt.1
    # ["190609_11_21_47.txt",     (11.9, 30.),     "sos" ,         "view"      ],  #Blade 7:Fatigue test, 41 hours at 4.48 Hz pt.1
    ["190613_16_08_36.txt",     (6, 20.),     "sos" ,         "save"      ],  #Blade 7:Fatigue test, 41 hours at 4.48 Hz pt.1
    
]





#Do not touch after this comment... all cases are setup up there..






from ver_datos import ver_datos
from get_datos import get_datos

def ver(fname,tmin,tmax,filter_type):
    ver_datos(fname,
        show_this_fig=True,
        save_this_fig=True,
        save_this_input=False,
        filter_type=filter_type,
        tmin=tmin,
        tmax=tmax,
        f1=f1,
        f2=f2)

def guardar(fname,tmin,tmax,filter_type,casename):
    ver_datos(fname,
        show_this_fig=False,
        save_this_fig=True,
        save_this_input=True,
        filter_type=filter_type,
        tmin=tmin,
        tmax=tmax,
        casename=casename,
        f1=f1,
        f2=f2)

with open("iga_inputs/00_start_and_end_times.txt", "w") as fid:
#              12345 12345 1234510 1234510789
    fid.write("#  t1    t2  case   fname\n")

    for i, case in enumerate(cases):
        fname = case[0]
        t1, t2 = case[1]
        filter_type = case[2]
        action = case[3]

        casename = "case{0:02.0f}".format(i)

        fid.write("{0:5} {1:5} {2:6} {3}\n".format(t1, t2, casename, fname))

        if action == "view" or action == "save" or action == "visa":
            get_datos(fname)
        if action == "view" or action == "visa":
            ver(fname,t1,t2,filter_type)
        if action == "save" or action == "visa":
            guardar(fname,t1,t2,filter_type,casename=casename)

