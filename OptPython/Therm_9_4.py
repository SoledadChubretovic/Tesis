def ThermalAnalysis(lambda_clay,b,ca,cb,cc,cd,dimz,n_cav):
    import numpy as np #para poder hacer matrices
    import sys
    from ansys.mapdl.core import launch_mapdl, Mapdl
    f = open(r"C:\Users\nchubretovic\OneDrive - Entel\Escritorio\Sole\OptPython\Error.txt", "a")

    c1a = ca[0:4] #m
    c2a = ca[4:8] #m
    c3a = ca[8:12] #m
    c4a = ca[12:16] #m
    c5a = ca[16:20] #m
    c6a = ca[20:24] #m
    c7a = ca[24:28] #m
    c8a = ca[28:32] #m
    c9a = ca[32:36] #m

    c1b = cb[0:4] #m
    c2b = cb[4:8] #m
    c3b = cb[8:12] #m
    c4b = cb[12:16] #m
    c5b = cb[16:20] #m
    c6b = cb[20:24] #m
    c7b = cb[24:28] #m
    c8b = cb[28:32] #m
    c9b = cb[32:36] #m

    c1c = cc[0:4] #m
    c2c = cc[4:8] #m
    c3c = cc[8:12] #m
    c4c = cc[12:16] #m
    c5c = cc[16:20] #m
    c6c = cc[20:24] #m
    c7c = cc[24:28] #m
    c8c = cc[28:32] #m
    c9c = cc[32:36] #m

    c1d = cd[0:4] #m
    c2d = cd[4:8] #m
    c3d = cd[8:12] #m
    c4d = cd[12:16] #m
    c5d = cd[16:20] #m
    c6d = cd[20:24] #m
    c7d = cd[24:28] #m
    c8d = cd[28:32] #m
    c9d = cd[32:36] #m

    try:
        mapdl = Mapdl(start_instance=False)
        # CALCULO TERMICO LADRILLO

        mapdl.finish()
        mapdl.clear()
        mapdl.prep7()

        # INPUTS
        matriz = np.array([c1a,c2a,c3a,c4a,c5a,c6a,c7a,c8a,c9a,c1b,c2b,c3b,c4b,c5b,c6b,c7b,c8b,c9b,c1c,c2c,c3c,c4c,c5c,c6c,c7c,c8c,c9c,c1d,c2d,c3d,c4d,c5d,c6d,c7d,c8d,c9d])

        ###### GEOMETRIA ######

        brick_perimeter = mapdl.blc5(b[0],b[1],b[2],b[3])
        for row in matriz:
            mapdl.blc5(row[0],row[1],row[2],row[3])

        cavities = mapdl.asba(brick_perimeter, 'all')
        mapdl.allsel()
        mapdl.vext(cavities,dz = dimz) #m

        ###### ASIGNACION DE ATRIBUTOS ###### 

        # Material (Conductividad termica, material isotrópico)
        Id_arcilla = 1
        mapdl.mp("KXX",Id_arcilla,lambda_clay)
        mapdl.mp("KYY",Id_arcilla,lambda_clay)
        mapdl.mp("KZZ",Id_arcilla,lambda_clay)
        
        # Elemento finito
        Id_ladrillo = 1
        mapdl.et(Id_ladrillo,"Solid70") #element reference p.297 "SOLID70"

        # Asignar todo a ladrillo
        mapdl.vatt(Id_arcilla,0,Id_ladrillo)

        ###### MALLADO ######
        mapdl.esize(0.006)
        mapdl.vsweep(1)
        
        ###### CONDICIONES DE BORDE DOF CONSTRAINS ######

        Thot = 30 #C
        Tcold = 18 #C

        mapdl.asel("S", vmin=2)
        mapdl.da("all", "temp", Tcold)
        mapdl.asel("S", vmin=4)
        mapdl.da("all", "temp", Thot)
        out = mapdl.allsel()

        ###### SOLVE ###### 
        mapdl.run("/SOLU")
        mapdl.solve()
        out = mapdl.finish()

        ###### POST-PROCESSING ###### 
        mapdl.post1()
        mapdl.set("last","last")

        ###### ITERACIONES CONVECCION ###### 

        iteracion = 0

        while iteracion < 2:
            # Calculo Tm
            # Necesito temperaturas en area 8,10,12,14,16,18,20,22
            a_conv = []
            for i in range(6,38):
                a_conv = np.append(a_conv,i)
            for i in range(39,n_cav*4+7):
                a_conv = np.append(a_conv,i)

            T_conv = np.zeros(len(a_conv))
            j = 0
            for i in  a_conv:
                mapdl.asel("S", vmin=i)
                mapdl.nsla()
                T = mapdl.post_processing.nodal_temperature() #temperaturas en cada nodo del area
                T_conv[j]= sum(T)/len(T) #temperatura promedio del area
                j = j+1

            Tm_conv = np.zeros(int(len(T_conv)/4)) #temperatura media para cada cavidad [c1,c2,c3,c4]
            j = 0
            for i in range(0,len(T_conv),4):
                Tm_conv[j] = (T_conv[i]+T_conv[i+1]+T_conv[i+2]+T_conv[i+3])/4
                j = j+1

            # Calculo parametro H (W/m2K)
            btz = 5.67e-8 #W/m2K4
            RAD = 0.83

            d = np.zeros(n_cav)
            for i in range(0,len(d)):
                d[i]= matriz[i,3] #m

            b1 = np.zeros(n_cav)
            for i in range(0,len(b1)):
                b1[i]= matriz[i,2] #m

            hr0 = 4*btz*(np.array(Tm_conv+273.15)**3)

            ha = np.zeros(n_cav)
            for i in range(0,len(ha)):
                ha[i] = max(0.025/d[i],1.25)

            hr = hr0/((1/RAD)+(1/RAD)-2+(2/(1+np.sqrt(1+(d**2)/(b1**2))-d/b1)))

            H = hr + ha #W/m2K

            # SURFACE LOAD CONVECTION

            mapdl.run("/SOLU")
            mapdl.sfcum("all","REPL") #reemplaza las cargas en areas

            j = 0
            for i in range(7,len(T_conv),4):
                mapdl.asel("S", vmin=i,vmax=i+3) #selecciona areas
                mapdl.nsla() #selecciona nodos asociados a las areas seleccionadas
                mapdl.sf("all","conv",H[j],Tm_conv[j])
                j = j+1

            out = mapdl.allsel()

            # SOLVE
            mapdl.solve()
            out = mapdl.finish()

            # POST-PROCESSING
            mapdl.post1()
            mapdl.set("last","last")
            #mapdl.post_processing.plot_nodal_temperature()
            
            iteracion = iteracion + 1


        ###### FLUJO DE CALOR ###### 

        # Cara fria
        mapdl.set("last","last")
        mapdl.asel("S", vmin=2)
        nodes = mapdl.nsla()

        min_nodenum_cold = int(mapdl.get("min_nodenum_cold","node","0","num","min")) #numero minimo de nodo seleccionado
        max_nodenum_cold = int(mapdl.get("max_nodenum_cold","node","0","num","max")) #numero maximo de nodo seleccionado
        nb_selected_nodes_cold = mapdl.mesh.n_node #numero de nodos seleccionados

        j = 0
        fcold = np.zeros(nb_selected_nodes_cold) 
        for i in range(min_nodenum_cold,max_nodenum_cold + 1):
            fcold[j] = mapdl.get("fcold","node",i,"tf","y")
            j = j + 1

        fcold = sum(abs(fcold))/len(fcold)

        # Cara caliente
        mapdl.set("last","last")
        mapdl.asel("S", vmin=4)
        nodes = mapdl.nsla()

        min_nodenum_hot = int(mapdl.get("min_nodenum_hot","node","0","num","min")) #numero minimo de nodo seleccionado
        max_nodenum_hot = int(mapdl.get("max_nodenum_hot","node","0","num","max")) #numero maximo de nodo seleccionado
        nb_selected_nodes_hot = mapdl.mesh.n_node #numero de nodos seleccionados

        j = 0
        fhot = np.zeros(nb_selected_nodes_hot) 
        for i in range(min_nodenum_hot,max_nodenum_hot + 1):
            fhot[j] = mapdl.get("fhot","node",i,"tf","y")
            j = j + 1

        fhot = sum(abs(fhot))/len(fhot)
        flux_area = (fcold+fhot)/2 #W/m2

        ###### CONDUCTIVIDAD EQUIVALENTE ###### 

        keq = flux_area*b[3]/(Thot-Tcold) #W/mK

        ###### TRANSMITANCIA DEL MURO ######

        RsiRse = 0.17 #m2K/W
        k_mort_pega = 0.23 #W/mK
        prop_ladrillo = 0.84
        prop_mort = 1 - prop_ladrillo
        C_muro = (prop_ladrillo*keq+prop_mort*k_mort_pega)/b[3] #W/m2K
        U_muro = 1/(RsiRse+1/C_muro) #W/m2K

        f.write("Ningun error. --- Umuro = "+ str(U_muro) + "\n")

    except:
        U_muro = 50
        f.write("Falla --- Umuro = "+ str(U_muro) + " Error:" )
        for i in range(len(sys.exc_info())):
            f.write(str(sys.exc_info()[i]) + "\n")
    
    f.close()
    
    return U_muro

