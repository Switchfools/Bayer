import numpy as np
import matplotlib.pyplot as plt
from numba import jit
from matplotlib.patches import Rectangle
from math import ceil
import matplotlib.animation as animation
from scipy.spatial.distance import pdist, squareform
from matplotlib import gridspec


class Persona():

    def __init__(self, x, y, xc, yc, grupo, trabajo, coordTrabajo,publicTransport, est=0):
        self.x = x
        self.y = y
        self.est = est
        # Sintomatico -5-4-3-2-1
        # Sano 0
        # Asintomatico 1 2 3 4 5
        self.grupo = grupo
        self.trabajo = trabajo
        self.coordTrabajo=coordTrabajo
        self.xc = xc
        self.yc = yc
        #0 if they have car and 1 if they have to use other type of transport
        self.publicTransport=publicTransport

    def moverseHacia(self,lugar):
        if(lugar=="Casa"):
            dist_x=self.x -self.xc
            dist_y=self.y -self.yc
            if(abs(dist_x)>abs(dist_y)):
                self.x-=np.sign(dist_x)
            else:
                self.y-=np.sign(dist_y)
        elif(lugar=="Trabajo"):
            dist_x=self.x -self.coordTrabajo[0]
            dist_y=self.y -self.coordTrabajo[1]
            if(abs(dist_x)>abs(dist_y)):
                self.x-=np.sign(dist_x)
            else:
                self.y-=np.sign(dist_y)


    def moverse(self):
        self.x=np.mod(self.x+np.random.randint(-1,2),self.coordTrabajo[2])
        self.y=np.mod(self.y+np.random.randint(-1,2),self.coordTrabajo[3])

    def cambiar(self,probSymptoms,probGetHeal,probDie):
        if (self.est> 0):
            if (np.random.rand()<probability(np.abs(self.est)*probSymptoms)):
                self.est=-1
            elif(np.random.rand()<probability(np.abs(self.est)*probGetHeal)):
                self.est=0
            else:
                self.est+=1
        elif (self.est < 0):
            if (np.random.rand()<probability(np.abs(self.est)*probGetHeal)):
                self.est=0
            elif(np.random.rand()<probability(np.abs(self.est)*probDie)):
                self.est=None
            else:
                self.est-=1
        return self.est

    def enfermarse(self,prob_getsick):
        if np.random.rand()<prob_getsick:
            self.est=1
        return(self.est)

    def EnCasa(self):
        en_casa=False
        if(self.x==self.xc and self.y==self.yc):
            en_casa=True
        return(en_casa)

    def EnTrabajo(self):
        en_trabajo=False
        if(self.x<=self.coordTrabajo[0]):
            en_trabajo=True
        return(en_trabajo)
class Mundo():

    def __init__(self, population,days, LaborSchedule,probs,size_x=1, size_y=1,asintomatics_x=0,asintomatics_y=0,sicks_x=0,sicks_y=0,day=0,hour=0,minute=0,Persons=[],Nhealt=0,Ninfected=0,Ndead=0,AbletoWork=0,MaximunPerDay=0):
        self.size_x = size_x
        self.size_y = size_y
        self.population = population
        ##the next parameters are statistical parameters
        self.Nhealt=population
        self.Ninfected=Ninfected
        self.Ndead=Ndead
        self.asintomatics_x= asintomatics_x
        self.asintomatics_y= asintomatics_y
        self.sicks_x = sicks_x
        self.sicks_y = sicks_y
        self.healthy_x = population
        self.healthy_y = population
        #Now we got Simulation parameters
        self.LaborSchedule=LaborSchedule
        self.probs=probs
        #Finally the number of days, that the simulation will go on
        self.days = days
        self.day=day
        self.hour=hour
        self.minute=minute
        #Personas en el mundo
        self.Persons=Persons
        #Personas Que no estan trabajando
        self.AbletoWork=AbletoWork
        self.MaximunPerDay=MaximunPerDay
    def CreateMap(self,NumberOfBuildigs, DimensionBuildings,estrategia):
        DimensionHouses=np.array([4,4])
        SizeBetweenBuildings= 4
        SpaceBetweenHouses=1
        DeadSpace=3*DimensionBuildings[0]
        self.size_y= NumberOfBuildigs*DimensionBuildings[1]+SizeBetweenBuildings*(NumberOfBuildigs+1)
        size_x_Houses= ceil((self.population/ceil(self.size_y/(DimensionHouses[1]+2*SpaceBetweenHouses)))*(DimensionHouses[0]+2*SpaceBetweenHouses))
        self.size_x= DimensionBuildings[0]+DeadSpace+size_x_Houses
        Map=np.zeros((self.size_x,self.size_y))
        Persons=[]
        #Now we fill the map
        ## First we put the buildings
        someX, someY=0,0
        currentAxis = plt.gca()
        currentAxis.set_axis_off()
        for i in range(NumberOfBuildigs):
            someX, someY = DimensionBuildings[0]/2, someY+SizeBetweenBuildings+(DimensionBuildings[1]/2)
            currentAxis.add_patch(Rectangle((someX-DimensionBuildings[0]/2, someY-DimensionBuildings[1]/2), DimensionBuildings[0], DimensionBuildings[1],fill=None))
            someX, someY = DimensionBuildings[0]/2, someY+DimensionBuildings[1]/2
        #Now we fill it with houses
        someX, someY = DimensionBuildings[0]+DeadSpace,SpaceBetweenHouses+DimensionHouses[1]/2
        counter=0
        for i in range(ceil(self.size_y/(DimensionHouses[1]+2*SpaceBetweenHouses))):
            someX,someY=DimensionBuildings[0]+DeadSpace,someY+DimensionHouses[1]/2
            for j in range(ceil((self.population/ceil(self.size_y/(DimensionHouses[1]+2*SpaceBetweenHouses))))):
                someX, someY = someX+DimensionHouses[0]/2+SpaceBetweenHouses, someY
                if(counter<500):
                    p=Persona(x=someX, y=someY, xc=someX, yc=someY, grupo=self.GroupsByestrategy(estrategia), trabajo=0,coordTrabajo=[np.random.randint(DimensionBuildings[0]),np.random.randint(self.size_y),DimensionBuildings[0],self.size_y], publicTransport=np.random.randint(2))
                    Persons.append(p)
                currentAxis.add_patch(Rectangle((someX-DimensionHouses[0]/2, someY-DimensionHouses[1]/2), DimensionHouses[0], DimensionHouses[1],fill=None))
                someX, someY = someX+DimensionHouses[0]/2, someY
                counter+=1
            someY+=DimensionHouses[1]/2+SpaceBetweenHouses
        plt.imshow(np.transpose(Map), cmap="cividis")
        plt.savefig("Map.png",bbox_inches='tight', pad_inches=0,transparent=True)
        self.Persons=Persons
        return(Persons)
    def GroupsByestrategy(self,estrategia):
        if(estrategia==0):
        #Progressive returns by groups
            group=np.random.choice(3, p=[60/500, (60+250)/500, (500-2*60-250)/500])
        elif(estrategia==1):
        #doing nothing at all, everyone comes back to work
            group=0
        elif(estrategia==2):
            group=np.random.choice(5)
        return(group)
    def Estrategia(self,estrategia,p):
        if(estrategia==0):
        #Progressive returns by groups
            if(self.day%7<5):
                if(p.grupo==0 and self.day>=0):
                    p.moverseHacia("Trabajo")
                elif(p.grupo==1 and self.day>=5):
                    p.moverseHacia("Trabajo")
                elif(p.grupo==2 and self.day>=10):
                    p.moverseHacia("Trabajo")
        elif(estrategia==1):
            #doing nothing at all, everyone comes back to work, but resting on weekends
            if((self.day%7)<5):
                p.moverseHacia("Trabajo")
        elif(estrategia==2):
            if(p.grupo==self.day%7):
                p.moverseHacia("Trabajo")




    def Step(self,estrategia):
        self.asintomatics_x= []
        self.asintomatics_y= []
        self.sicks_x = []
        self.sicks_y = []
        self.healthy_x = []
        self.healthy_y = []
        if(self.hour<self.LaborSchedule[0] or self.hour>self.LaborSchedule[1]):
            #people at either at home or moving
            self.AbletoWork=0
            for p in self.Persons:
                if(p.est==0):
                    self.healthy_x.append(p.x)
                    self.healthy_y.append(p.y)
                    if(p.EnCasa()):
                        est=p.enfermarse(probability(self.probs[4]))
                        self.Ninfected+=est
                        self.Nhealt-=est
                    else:
                        est=p.enfermarse(probability(self.probs[5]*p.publicTransport))
                        p.moverseHacia("Casa")
                        self.Ninfected+=est
                        self.Nhealt-=est
                elif(p.est>0):
                    self.asintomatics_x.append(p.x)
                    self.asintomatics_y.append(p.y)
                    state=p.cambiar(probs[0],probs[1],probs[2])
                    if(state==None):
                        self.Persons.remove(p)
                        self.Ninfected-=1
                        self.Ndead+=1
                    elif(state==0):
                        self.Ninfected-=1
                        self.Nhealt+=1

                    if not (p.EnCasa()):
                        p.moverseHacia("Casa")

                else:
                    self.sicks_x.append(p.x)
                    self.sicks_y.append(p.y)
                    state=p.cambiar(probs[0],probs[1],probs[2])
                    if(state==None):
                        self.Persons.remove(p)
                        self.Ninfected-=1
                        self.Ndead+=1
                    elif(state==0):
                        self.Ninfected-=1
                        self.Nhealt+=1
                    if not (p.EnCasa()):
                        p.moverseHacia("Casa")
        else:
            positions=[]
            infecteds=[]
            for p in self.Persons:
                positions.append([p.x,p.y])
                if(p.est!=0):
                    infecteds.append(1)
                else:
                    infecteds.append(0)
            D = squareform(pdist(positions))
            i=0
            self.AbletoWork=0
            for p in self.Persons:
                if(p.est==0):
                    self.healthy_x.append(p.x)
                    self.healthy_y.append(p.y)
                    if(p.EnCasa()):
                        est=p.enfermarse(probability(self.probs[4]))
                        self.Estrategia(estrategia,p)
                        self.Ninfected+=est
                        self.Nhealt-=est
                    elif(p.EnTrabajo()):
                        p.moverse()
                        cercaInfectado=np.less(D[i,infecteds],1.8)
                        if(np.sum(cercaInfectado)>0):
                            #hace falta calcular las distancias
                            est=p.enfermarse(self.probs[3])
                            self.Ninfected+=est
                            self.Nhealt-=est
                    else:
                        est=p.enfermarse(probability(self.probs[5]*p.publicTransport))
                        p.moverseHacia("Trabajo")
                        self.Ninfected+=est
                        self.Nhealt-=est
                elif(p.est>0):
                    self.asintomatics_x.append(p.x)
                    self.asintomatics_y.append(p.y)
                    state=p.cambiar(probs[0],probs[1],probs[2])
                    if(state==None):
                        self.Persons.remove(p)
                        self.Ninfected-=1
                        self.Ndead+=1
                    elif(state==0):
                        self.Ninfected-=1
                        self.Nhealt+=1
                    if(p.EnTrabajo()):
                        p.moverse()
                    elif (p.EnCasa()):
                        self.Estrategia(estrategia,p)
                    else:
                        p.moverseHacia("Trabajo")

                else:
                    self.sicks_x.append(p.x)
                    self.sicks_y.append(p.y)
                    state=p.cambiar(probs[0],probs[1],probs[2])
                    if(state==None):
                        self.Persons.remove(p)
                        self.Ninfected-=1
                        self.Ndead+=1
                    elif(state==0):
                        self.Ninfected-=1
                        self.Nhealt+=1
                    if not (p.EnCasa()):
                        p.moverseHacia("Casa")
                if(p.EnTrabajo()):
                    self.AbletoWork+=1
                i+=1
        self.minute+=1
        if(self.AbletoWork>self.MaximunPerDay):
            self.MaximunPerDay=self.AbletoWork
        if(self.minute==60):
            self.minute=0
            self.hour+=1
            #return(self.Nhealt,self.Ninfected,self.Ndead)
            if(self.hour==24):
                self.hour=0
                self.day+=1
                return(self.Nhealt,self.Ninfected,self.Ndead,self.MaximunPerDay)
                self.MaximunPerDay=0
                for p in self.Persons:
                    if(p.est>0):
                        p.est+=1
                    elif(p.est<0):
                        p.est-=1
            #people at either at wrok, moving, or in their own houses depending on the group


def probability(rate):
    #the rate is given in days, we computed in minutes
    rate_min=rate/(24*60)
    return(1-np.exp(-rate_min))



def init():
    """initialize animation"""
    global MiMundo
    Healthy.set_data([], [])
    Assimptomatics.set_data([], [])
    Simptomatics.set_data([], [])
    Sanos.set_data([], [])
    Infectados.set_data([], [])
    Muertos.set_data([], [])
    Trabajan.set_data([],[])
    return Healthy,Assimptomatics,Simptomatics,Sanos,Infectados,Muertos,Trabajan
def animate(i):
    """perform animation step"""
    global MiMundo,Persons, ax, fig,DataCurve
    a=MiMundo.Step(estrategia=1)
    if(a!=None):
        S.append(a[0])
        I.append(a[1])
        M.append(a[2])
        T.append(a[3])
    # update pieces of the animation
    Healthy.set_data(MiMundo.healthy_x, MiMundo.healthy_y)
    Assimptomatics.set_data(MiMundo.asintomatics_x, MiMundo.asintomatics_y)
    Simptomatics.set_data(MiMundo.sicks_x, MiMundo.sicks_y)
    Sanos.set_data(range(len(S)), S)
    Infectados.set_data(range(len(I)),I)
    Muertos.set_data(range(len(M)), M)
    Trabajan.set_data(range(len(M)),T)
    return Healthy,Assimptomatics,Simptomatics,Sanos,Infectados,Muertos,Trabajan


if __name__ == "__main__":
    #First We Define the Parameters
    LaborSchedule=(9,16)
    population=500
    days=30
    #this is an array for the rates related to the probabilities [probSymptoms,probGetHeal,probDie,prob_getsick, prob_getsick_at_home,prob_getsick_in_public_transport]
    LongTermSicks=(0.0000001*population)
    probs=[(0.4*LongTermSicks)/days,(0.2514*LongTermSicks)/days,(0.000621*LongTermSicks)/days,0.3,(0.05*LongTermSicks)/days, (0.6*LongTermSicks)/days]
    print(probability(np.array(probs)))
    #MiMundo=Mundo(population=population,days=days,LaborSchedule=LaborSchedule,probs=probs)
    #We Strat the simulation
    #MiMundo.SimulateWithDraw(Persons,LaborSchedule,probs)
    '''
    fig = plt.figure()
    spec = gridspec.GridSpec(ncols=1, nrows=2,height_ratios=[5,2])
    ax = fig.add_subplot(spec[0])
    Persons=MiMundo.CreateMap(3,(10,40),estrategia=1)
    S=[]
    I=[]
    M=[]
    T=[]
    Healthy, = ax.plot([], [], 'go', ms=3)
    Assimptomatics, = ax.plot([], [], 'bo', ms=3)
    Simptomatics, = ax.plot([], [], 'ro', ms=3)
    ax2 = fig.add_subplot(spec[1])
    ax2.set_xlim([0,(days*24)*60])
    ax2.set_ylim([0, population])

    Sanos, = ax2.plot([], [],lw=3,c="teal")
    Infectados, = ax2.plot([], [],lw=3,c="orange")
    Muertos, = ax2.plot([], [],lw=3,c="red")
    Trabajan, = ax2.plot([], [],lw=1,c="purple", alpha=0.5)
    #ax.add_patch(pa)
    #ax1 = fig.add_subplot(1,1,1)
    ani = animation.FuncAnimation(fig, animate, frames=60,interval=1, blit=True, init_func=init)
    plt.show()
    print(M[-1])
    '''
    Samples=3
    s=np.zeros(30)
    i=np.zeros(30)
    m=np.zeros(30)
    t=np.zeros(30)
    for experiment in range(Samples):
        S=[]
        I=[]
        M=[]
        T=[]
        MiMundo=Mundo(population=population,days=days,LaborSchedule=LaborSchedule,probs=probs)
        Persons=MiMundo.CreateMap(3,(10,40),estrategia=1)
        while(MiMundo.day<days):
            a=MiMundo.Step(estrategia=1)
            if(a!=None):
                S.append(a[0])
                I.append(a[1])
                M.append(a[2])
                T.append(a[3])
        s=s+np.array(S)
        i=i+np.array(I)
        m=m+np.array(M)
        t=t+np.array(T)
        print("experimento #:",experiment)

    print("finish")
    plt.figure()
    dias=range(30)
    plt.plot(dias,s/Samples,c="green", label='Sanos')
    plt.plot(dias,i/Samples,c="orange", label='Infectados')
    plt.plot(dias,m/Samples,c="red", label='muertos')
    plt.plot(dias,t/Samples,c="teal", label='Trabajando')
    plt.legend()
    plt.grid()
    plt.title("Promedio Simulaciones para estrategia 1")
    plt.savefig("Experimento_Estrategia_1.pdf")
    plt.show()
