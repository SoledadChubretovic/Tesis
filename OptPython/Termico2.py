def ThermalAnalysis(tmin,k,b,c1,c2,c3,c4,dimz):
        
    # Unidades: mm,W,C

    import numpy as np #para poder hacer matrices
    import sys
    from ansys.mapdl.core import launch_mapdl, Mapdl
    f = open(r"C:\Users\nchubretovic\OneDrive - Entel\Escritorio\Sole\OptPython\Error.txt", "a")

    if (c1[2]+c2[2]+3*tmin) - b[2] > 0 or (c3[2]+c4[2]+3*tmin) - b[2] > 0 or (c1[3]+c3[3]+3*tmin) - b[3] > 0 or (c2[3]+c4[3]+3*tmin) - b[3] > 0:
        U_muro =  50
        f.write("Falla por cavidades fuera. --- Umuro = "+ str(U_muro) + "\n" )
    elif c1[3]+c4[3]+3*tmin > b[3] and c1[2]+c4[2]+3*tmin > b[2]:
        U_muro =  50
        f.write("Falla por traslapo 1. --- Umuro = "+ str(U_muro) + "\n" )
    elif c2[3]+c3[3]+3*tmin > b[3] and c2[2]+c3[2]+3*tmin > b[2]:
        U_muro =  50
        f.write("Falla por traslapo 2. --- Umuro = "+ str(U_muro) + "\n" )
    else:
        try:
            mapdl = Mapdl(start_instance=False)
            # CALCULO TERMICO LADRILLO

            mapdl.finish()
            mapdl.clear()
            mapdl.prep7()

            # INPUTS
            matriz = np.array([c1,c2,c3,c4])

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
            mapdl.mp("KXX",Id_arcilla,k)
            mapdl.mp("KYY",Id_arcilla,k)
            mapdl.mp("KZZ",Id_arcilla,k)
            
            # Elemento finito
            Id_ladrillo = 1
            mapdl.et(Id_ladrillo,"Solid70") #element reference p.297 "SOLID70"

            # Asignar todo a ladrillo
            mapdl.vatt(Id_arcilla,0,Id_ladrillo)

            ###### MALLADO ######
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
            #mapdl.post_processing.plot_nodal_temperature()

            ###### ITERACIONES CONVECCION ###### 

            iteracion = 0

            while iteracion < 7:
                # Calculo Tm
                # Necesito temperaturas en area 8,10,12,14,16,18,20,22
                a_conv = [8,10,12,14,16,18,20,22] #que coincida con el orden de la geometria
                T_conv = np.zeros(len(a_conv))
                j = 0

                for i in a_conv:
                    mapdl.asel("S", vmin=i)
                    mapdl.nsla()
                    T = mapdl.post_processing.nodal_temperature() #temperaturas en cada nodo del area
                    T_conv[j]= sum(T)/len(T) #temperatura promedio del area
                    j = j+1

                Tm_conv = np.zeros(int(len(T_conv)/2)) #temperatura media para cada cavidad [c1,c2,c3,c4]
                j = 0

                for i in range(0,len(T_conv),2):
                    Tm_conv[j] = (T_conv[i]+T_conv[i+1])/2
                    j = j+1

                # Calculo parametro H (W/m2K)
                btz = 5.67e-8 #W/m2K4
                RAD = 0.83

                d = np.zeros(4)
                for i in range(0,len(d)):
                    d[i]= matriz[i,3] #m

                b1 = np.zeros(4)
                for i in range(0,len(b1)):
                    b1[i]= matriz[i,2] #m

                hr0 = 4*btz*(np.array(Tm_conv+273.15)**3)

                ha = np.zeros(4)
                for i in range(0,len(ha)):
                    ha[i] = max(0.025/d[i],1.25)

                hr = hr0/((1/RAD)+(1/RAD)-2+(2/(1+np.sqrt(1+(d**2)/(b1**2))-d/b1)))

                H = hr + ha #W/m2K

                # SURFACE LOAD CONVECTION

                mapdl.run("/SOLU")
                mapdl.sfcum("all","REPL") #reemplaza las cargas en areas

                #cavidad 1
                mapdl.asel("S", vmin=7,vmax=10) #selecciona areas
                mapdl.nsla() #selecciona nodos asociados a las areas seleccionadas
                mapdl.sf("all","conv",H[0],Tm_conv[0])

                #cavidad 2
                mapdl.asel("S", vmin=11,vmax=14) #selecciona areas
                mapdl.nsla() #selecciona nodos asociados a las areas seleccionadas
                mapdl.sf("all","conv",H[1],Tm_conv[1])

                #cavidad 3
                mapdl.asel("S", vmin=15,vmax=18) #selecciona areas
                mapdl.nsla() #selecciona nodos asociados a las areas seleccionadas
                mapdl.sf("all","conv",H[2],Tm_conv[2])

                #cavidad 4
                mapdl.asel("S", vmin=19,vmax=22) #selecciona areas
                mapdl.nsla() #selecciona nodos asociados a las areas seleccionadas
                mapdl.sf("all","conv",H[3],Tm_conv[3])

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
            prop_mort = 1-prop_ladrillo
            C_muro = (prop_ladrillo*keq+prop_mort*k_mort_pega)/b[3] #W/m2K
            U_muro = 1/(RsiRse+1/C_muro) #W/m2K

            f.write("Ningun error. --- Umuro = "+ str(U_muro) + "\n")

        except:
            U_muro = 50
            f.write("Falla por traslapo 3. --- Umuro = "+ str(U_muro) + "\n" )
            for i in range(len(sys.exc_info())):
                f.write(str(sys.exc_info()[i]))
        
        f.close()
    
    return U_muro

