import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from scipy import interpolate        
from numpy import linspace
import numpy as np
import matplotlib.tri as tri
plt.rcParams['text.usetex'] = True


def plot_tricontourc( h_n, nodeid, centerc, vertexn):
       
    plt.figure(figsize=(8,4))
    x = centerc[:,0]
    y = centerc[:,1]
    #Hinterp  = interpolate.NearestNDInterpolator((x, y) , h_n)
    XI  = vertexn[:,0]
    YI  = vertexn[:,1] 
    #hI  = Hinterp(XI,YI)
    #plt.tricontour(XI, YI, h_n, 40)
    plt.tricontourf(x, y, h_n,60)#, cmap=plt.cm.jet)
    plt.title(r' $h+b$ at $t = 1$ ', fontsize=20)
    plt.xlim([0, 2])
    plt.ylim([0, 1])
    plt.savefig("hbc1.png")
    plt.show()
def plot_tricontourn( h_n, nodeid, centerc, vertexn):
       
    plt.figure(figsize=(8,4))
    x = centerc[:,0]
    y = centerc[:,1]
    #Hinterp  = interpolate.NearestNDInterpolator((x, y) , h_n)
    XI  = vertexn[:,0]
    YI  = vertexn[:,1] 
    #hI  = Hinterp(XI,YI)
    #plt.tricontour(XI, YI, h_n, 40)
    plt.tricontourf(XI, YI, h_n, cmap=plt.cm.jet)
    plt.title(r' $h+b$ at $t = 0.3$ ', fontsize=20)
    plt.xlim([0, 2])
    plt.ylim([0, 1])
    plt.savefig("grhbn1.png")
    plt.show()


def plot_tricontourcj( h_n, nodeid, centerc, vertexn):
       
    plt.figure(figsize=(8,4))
    x = centerc[:,0]
    y = centerc[:,1]
    #Hinterp  = interpolate.NearestNDInterpolator((x, y) , h_n)
    XI  = vertexn[:,0]
    YI  = vertexn[:,1] 
    #hI  = Hinterp(XI,YI)
    #plt.tricontour(XI, YI, h_n, 40)
    plt.tricontourf(x, y, h_n, 100, cmap=plt.cm.jet)
    plt.title(r' $h+b$ at $t = 0.3$ ', fontsize=20)
    plt.xlim([0, 2])
    plt.ylim([0, 1])
    plt.savefig("hb.png")
    plt.show()
def plot_tricontournj( h_n, nodeid, centerc, vertexn):
       
    plt.figure(figsize=(8,4))
    x = centerc[:,0]
    y = centerc[:,1]
    #Hinterp  = interpolate.NearestNDInterpolator((x, y) , h_n)
    XI  = vertexn[:,0]
    YI  = vertexn[:,1] 
    #hI  = Hinterp(XI,YI)
    #plt.tricontour(XI, YI, h_n, 40)
    plt.tricontourf(XI, YI, h_n, 100, cmap=plt.cm.jet)
    plt.title(r' $h+b$ at $t = 1$ ', fontsize=20)
    plt.xlim([0, 2])
    plt.ylim([0, 1])
    plt.savefig("jhbn1.png")
    plt.show()
               
         
         





def plot_tricontour( h_n, nodeid, centerc, vertexn, a):
       
    plt.figure()#figsize=(8,4))
    x = centerc[:,0]
    y = centerc[:,1]
   # Hinterp  = interpolate.NearestNDInterpolator((x, y) , h_n)
    XI  = vertexn[:,0]
    YI  = vertexn[:,1] 
    #hI  = Hinterp(XI,YI)
    #plt.tricontour(XI, YI, h_n, 40)
    if a == 1 :
         plt.tricontour(XI, YI, h_n, 40, cmap=plt.cm.jet)
         plt.title(r' $h$ ')
         plt.savefig("hn.png")
         plt.show()
         
    elif a == 2 :
         plt.tricontour(XI, YI, h_n, 40, cmap=plt.cm.jet)
         plt.title(r' $\vert \mathbf{u} \vert$ ') 
         plt.savefig("un.png")   
         plt.show()
         

    elif a == 3 :
         plt.tricontour(XI, YI, h_n, 40, cmap=plt.cm.jet)
         plt.title(r' $\vert \mathbf{B} \vert$  ')
         plt.savefig("Bn.png")
         plt.show()
         


    
def plot_tricontourdvv( h_n, nodeid, centerc, vertexn, a):
       
    plt.figure()#figsize=(8,4))
    x = centerc[:,0]
    y = centerc[:,1]
   # Hinterp  = interpolate.NearestNDInterpolator((x, y) , h_n)
    XI  = vertexn[:,0]
    YI  = vertexn[:,1] 
    #hI  = Hinterp(XI,YI)
    #plt.tricontour(XI, YI, h_n, 40)
    if a == 1 :
         plt.tricontour(XI, YI, h_n, 40, cmap=plt.cm.jet)
         plt.title(r' $h$ ')
         plt.savefig("hn.png")
         plt.show()
         
    elif a == 2 :
         plt.tricontour(XI, YI, h_n, 40, cmap=plt.cm.jet)
         plt.title(r' $\vert \mathbf{u} \vert$ ') 
         plt.savefig("un.png")   
         plt.show()
         

    elif a == 3 :
         plt.tricontour(XI, YI, h_n, 40, cmap=plt.cm.jet)
         plt.title(r' $\vert \mathbf{B} \vert$  ')
         plt.savefig("Bn.png")
         plt.show()
         
 
    
def plot_tripcolor( h_n, nodeid, centerc, vertexn):
       
    plt.figure(figsize=(7,7))
    x = centerc[:,0]
    y = centerc[:,1]
    fig, ax = plt.subplots()
    #ax.set_aspect('equal')
    TRI = tri.Triangulation(x, y)
    TP = ax.tripcolor(TRI, h_n, cmap=plt.cm.jet)
    fig.colorbar(TP)     
    plt.savefig("tricontour.png") 
    plt.show()
    plt.savefig("tricontour.png")   
    

#
#TP = ax.tripcolor(TRI, B2)
#TP = ax.tricontourf(TRI, B2)
#ax.tricontour(TRI, B2)
#fig.colorbar(TP)
def plot_quiver(h_c, hu_c, hv_c, h_n, hu_n, nodeid, centerc, vertexn, time):
       
    plt.figure(figsize=(7,7))
    x = centerc[:,0]
    y = centerc[:,1]
    plt.quiver(x, y, hu_c/h_c, hv_c/h_c)
    plt.show()
    plt.savefig("quiver.png")
       
    fig = plt.figure(figsize=(9,6))
    ax = fig.gca(projection='3d')  
    collec = ax.plot_trisurf(vertexn[:,0], vertexn[:,1], h_n, linewidth=2, antialiased=True,
                            cmap=cm.jet , triangles=nodeid)
    collec.autoscale()
    ax.view_init(20, -45)
    ax.grid(b=False)
    plt.show()
    plt.savefig("trisurf.png")
       
    XI = linspace(min(vertexn[:,0]), max(vertexn[:,0]), 100)
    YI = 125*np.ones(len(XI))
    Hinterp = interpolate.NearestNDInterpolator((vertexn[:,0], vertexn[:,1]), h_n)
    Uinterp = interpolate.NearestNDInterpolator((vertexn[:,0], vertexn[:,1]), hu_n/h_n)
       
    hI = Hinterp(XI,YI)
    plt.figure(figsize=(7,7))
    plt.plot(XI, hI)
       
    uI = Uinterp(XI,YI)
    plt.figure(figsize=(7,7))
    plt.plot(XI, uI)    
    
def plot_cross_section(h_n, hu_n, nodeid, centerc, vertexn):
       
    plt.figure(figsize=(7,7))
    x = centerc[:,0]
    y = centerc[:,1]

       
    XI = linspace(min(vertexn[:,0]), max(vertexn[:,0]), 10000)
    YI = 125*np.ones(len(XI))
    Hinterp = interpolate.NearestNDInterpolator((vertexn[:,0], vertexn[:,1]), h_n)
    Uinterp = interpolate.NearestNDInterpolator((vertexn[:,0], vertexn[:,1]), hu_n/h_n)
       
    hI = Hinterp(XI,YI)
    plt.figure(figsize=(7,7))
    plt.plot(XI, hI)
    plt.savefig("hn.png") 
       
    uI = Uinterp(XI,YI)
    plt.figure(figsize=(7,7))
    plt.plot(XI, uI)   
        
    plt.savefig("un.png")   
    plt.show() 
    
    

def plot_cross_sectionF(h_c, hu_c, hv_c, h_n, hu_n, nodeid, centerc, vertexn, time):
       
    plt.figure(figsize=(7,7))
    x = centerc[:,0]
    y = centerc[:,1]

       
    XI = linspace(min(vertexn[:,0]), max(vertexn[:,0]), 100)
    YI = 125*np.ones(len(XI))
    Hinterp = interpolate.NearestNDInterpolator((vertexn[:,0], vertexn[:,1]), h_n)
    Uinterp = interpolate.NearestNDInterpolator((vertexn[:,0], vertexn[:,1]), hu_n/h_n)
       
    hI = Hinterp(XI,YI)
    plt.figure(figsize=(7,7))
    plt.plot(XI, hI)
       
    uI = Uinterp(XI,YI)
    plt.figure(figsize=(7,7))
    plt.plot(XI, uI)    
        
    
    
