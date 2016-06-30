# -*- coding: utf-8 -*-
"""
Created on Mon Jun 27 19:32:23 2016

@author: YinCheang
"""

import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as sp

class gasdy :
    
    def Mayerfun(self,M):
        
        f = np.sqrt((2.4/0.4))*np.arctan(np.sqrt(0.4/2.4)*np.sqrt(M**2 -1)) - np.arctan(np.sqrt(M**2 -1))        
        
        return f
    
    def cal_Mayerfun(self,v):

        Guess = np.array([1.5])        
        
        f = lambda M : v-(np.sqrt((2.4/0.4))*np.arctan(np.sqrt(0.4/2.4)*np.sqrt(M**2 -1))- np.arctan(np.sqrt(M**2 -1)))
                        
        Mans = sp.fsolve(f,Guess)   
        
        return Mans
        
    def cal_MaAng(self,M):
    
        MaAng = np.arcsin(1./M)
        
        return MaAng
        
        
class Nozzle_Design(gasdy):
    
    Ts    = 26.5*np.pi/180.
    Ts_in = Ts/2
    
    Nwav  = 10
    slop  = Ts_in/Nwav
    ntime = 0.0

    wallx = []
    wally = []
    wallv = []
    walla = []

    pointx = []
    pointy = []
    pointv = [] 
    pointa = []
    
    sysx = []
    sysy = []    
    
    
    def cal_angle(self,M1,M2,i):        

#        ans = ( self.cal_MaAng(M1)-i*(self.slop) + (self.cal_MaAng(M2)-(i+1)*self.slop) ) /2
        ans = ( self.cal_MaAng(M1) + (self.cal_MaAng(M2)-self.slop) ) /2
        print
        print
        print "............cal_angel............"
        print M1 , M2
        print (180/np.pi)*self.cal_MaAng(M1)#-(180/np.pi)*i*(self.slop)
        print (180/np.pi)*ans , i 
        print (180/np.pi)*self.cal_MaAng(M2)#-(180/np.pi)*(i+1)*self.slop
        print 
        print "............cal_angel............"
        print
        print        
        
        return ans
    
    def inishape(self):
        
        self.wallx.append(0.0)
        self.wally.append(1.0)
        self.wallv.append(0)
        self.walla.append(0.0)
               
#=============================================================================
               
        self.wallx.append(0.5)
        self.wally.append(1.0)

        M1  = self.cal_Mayerfun(0.0)
        M2  = self.cal_Mayerfun(self.slop)             
        ans = self.cal_angle(M1,M2,0)
        
        
        self.walla.append(ans)

        
        M   = 1./np.sin(ans)
        ans = self.Mayerfun(M) # this is MayerFun V
        
        print ans 
        

        self.wallv.append(ans)#ans)
        self.pointx.append(0.5)
        self.pointy.append(1.0)
        self.pointv.append(ans)
        self.pointa.append(ans)
        
        
        intervalx = 0.3/(self.Nwav-1)    
        
        
        for i in range(self.Nwav-1):
            
#============================================================================
            pt = len(self.wallx)-1            
            
            x = self.wallx[pt] + intervalx
            y = np.tan( (i+1)*self.slop )*(x-self.wallx[pt]) + self.wally[pt]    
            self.wallx.append(x)
            self.wally.append(y)
#============================================================================
            
            
            v = (i+1)*self.slop   
            
            print            
            print
            print "Mayer v "
            print v*180/np.pi , (v+ self.slop)*180/np.pi 
            print
            print


            M1  =self.cal_Mayerfun(v)
            M2  =self.cal_Mayerfun(v + self.slop)
            
            print "M1 , M2"
            print M1 , M2 
            print
    
            ans = self.cal_angle(M1,M2,i+1)
            self.walla.append(ans)
            
            
#            print ans*180/np.pi
#            self.wallv.append(v)            
#============================================================================            
            ans = 1./np.sin(ans)
            print "M "
            print ans 
            
            ans = self.Mayerfun(ans)
            print
            print "Mayer v"
            print ans*180/np.pi
            print 
            
            self.wallv.append(ans)              
            
            self.pointx.append(x)
            self.pointy.append(y)
            self.pointv.append(ans)
            self.pointa.append(ans)
            
        for i in range(self.Nwav):
            self.ntime = self.ntime + (i+1)
            


#        pause = raw_input("wait")
        
#        plt.plot(self.wallx,self.wally,'b*')
#        plt.plot(self.pointx,self.pointy)
#        plt.show()
        
        return 0 
    
    def ini_caculation(self,x1,y1,x2,y2,m1,m2):
        
        print ".........m1...m2..........."
        print x1 , y1
        print x2 , y2
        print m1*180/np.pi , m2*180/np.pi
        print ".........m1...m2..........."    
        
        m1 = -np.tan(m1)
        m2 =  np.tan(m2)
        
        a = np.array([ [1.0,-m1]
                      ,[1.0,-m2] ] )
        b = np.array([ y1-m1*x1
                      ,y2-m2*x2] )
        
        a = np.linalg.inv(a)
        ans = np.dot(a,b)

        y   = ans[0]
        x   = ans[1]
            
        return x,y 
        
    def sys_caculation(self,m,x1,y1):
        
        print m*180/np.pi
        m = np.tan(m)
        print m 
        
        y = 0.0
        x = (-1./m)*(y-y1) + x1
        
        self.sysx.append(x)
        self.sysy.append(0.0)
        self.pointx.append(x)
        self.pointy.append(y)
        
        print x , y
        
        return 0
    def wall_caculation(self):

        cco = 0 
        pts  = 2*self.Nwav-1
        
        for i in range(int(self.Nwav)) :
#            print i , pts+cco , cco
            pt = len(self.wallx)-1
                
            x1  = self.wallx[pt]
            y1  = self.wally[pt]
                
            x2  = self.pointx[pts+cco]
            y2  = self.pointy[pts+cco]
                
            m1  = self.Ts_in - (i)*self.slop
                
            M  = self.cal_Mayerfun(self.pointv[pts+cco])
            m2 = self.cal_MaAng(M)  + self.pointa[pts+cco]           
            
            pts  = pts +cco
            cco = cco-1 
                
            x,y = self.ini_caculation(x1,y1,x2,y2,-1*m1,m2)
                
            self.pointx.append(x)
            self.pointy.append(y)   
            self.wallx.append(x)
            self.wally.append(y)
            
            if i == 0 :
                cco = self.Nwav-1
                
#            plt.plot(self.wallx,self.wally,'b*')
#            plt.plot(self.wallx,self.wally,'b')             
#            plt.plot(self.pointx,self.pointy,'r*')
#            plt.show()
            
#            pause = raw_input("WAIT")
            
        return 0    
            
    def caculation(self):
        
        cco = 0
        pcc = 0
##        swt = 3

        
        v = lambda v1,v2,s1,s2 : 0.5*(v1+v2)+0.5*(s1-s2)
        a = lambda v1,v2,s1,s2 : 0.5*(v1-v2)+0.5*(s1+s2)
        
        for i in range(int(self.Nwav),int(self.ntime+self.Nwav)):
            
            if cco == 0 :
                x1 = self.pointx[i-self.Nwav+pcc]
                y1 = self.pointy[i-self.Nwav+pcc]
                
                M  = self.cal_Mayerfun(self.pointv[i-self.Nwav+pcc ])
                print "==========M==========="
                print M
                print
                m  = self.cal_MaAng(M)- self.pointa[i-self.Nwav+pcc ]              
                
##                if swt > 0 :
##                    m = self.walla[i-3]
                
#                print i , i-self.Nwav+pcc   , m
                self.sys_caculation(m,x1,y1)
                
#==============================================================================                
                Q     = self.pointv[i-self.Nwav+pcc ] + self.pointa[i-self.Nwav+pcc ]       
                swall = 0.0
                vwall = Q-swall            
                
                print vwall , Q 
                temp = self.cal_Mayerfun(vwall)
                print temp , self.cal_MaAng(temp)*180/np.pi
                
#                print self.pointv[i-self.Nwav+pcc ]*180/np.pi , self.pointa[i-self.Nwav+pcc ]*180/np.pi
#==============================================================================          
                
                self.pointv.append(vwall)
                self.pointa.append(0.0)
                
                pcc = pcc + 1
                cco = self.Nwav - pcc         
                
#                plt.plot(self.wallx,self.wally,'b*')
#                plt.plot(self.wallx,self.wally,'b')             
#                plt.plot(self.pointx,self.pointy,'r*')
#                plt.show()
#                pause = raw_input("waittt")
                
                
            else:
                pdd = pcc - 1
                

                
                x1 = self.pointx[i-self.Nwav+pdd]
                y1 = self.pointy[i-self.Nwav+pdd] 
                x2 = self.pointx[i-1]
                y2 = self.pointy[i-1]

                M  = self.cal_Mayerfun(self.pointv[i-self.Nwav+pdd])
                print " M "
                print M
                m1  = self.cal_MaAng(M)-self.pointa[i-self.Nwav+pdd]   

##                if swt > 0 :
##                    m1 = self.walla[i-3]                           
                                                        
                M  = self.cal_Mayerfun(self.pointv[i-1])
                print M
                m2 = self.cal_MaAng(M)+self.pointa[i-1]

                print i, i-self.Nwav+pdd, i-1         
                
                x,y = self.ini_caculation(x1,y1,x2,y2,m1,m2)
                print x, y 
                
                tv = v(self.pointv[int(i-self.Nwav+pdd)],self.pointv[i-1],self.pointa[int(i-self.Nwav+pdd)],self.pointa[i-1])
                ta = a(self.pointv[int(i-self.Nwav+pdd)],self.pointv[i-1],self.pointa[int(i-self.Nwav+pdd)],self.pointa[i-1])
                
                self.pointx.append(x)
                self.pointy.append(y)                
                self.pointv.append(tv)
                self.pointa.append(ta)
                cco = cco - 1

                
#                plt.plot(self.wallx,self.wally,'b*')
#                plt.plot(self.wallx,self.wally,'b')             
#                plt.plot(self.pointx,self.pointy,'r*')
#                plt.show()
#                pause = raw_input("wait")  
                
##            swt= swt -1
                
        
        return 0
    
    def display(self):
        
        ax = []
        ay = []
        
        j = 0
        for i in range(self.Nwav):
            
            ax.append(self.wallx[i+1])
            ax.append(self.sysx[i])
            ay.append(self.wally[i+1])
            ay.append(self.sysy[i])
            plt.plot(ax[j:j+2],ay[j:j+2])
            j = j+2
            
        
        num = self.Nwav
        for i in range(num,-1,-1):
            plt.plot(self.pointx[num:num+i],self.pointy[num:num+i])
            num = num+i

        bx = []
        by = []
        
        j = 0
        pt = 2*self.Nwav-1
        ps = self.Nwav+1
        cco = 0
        for i in range(self.Nwav):
            print pt+cco
            bx.append(self.pointx[pt+cco])
            bx.append(self.wallx[ps+i])
            by.append(self.pointy[pt+cco])
            by.append(self.wally[ps+i])
                        

            plt.plot(bx[j:j+2],by[j:j+2])   
            j = j+2
            
            pt = pt+cco
            
            if i == 0:
                cco = self.Nwav
            
            cco = cco-1

                
        plt.plot(self.wallx,self.wally,'b*')
        plt.plot(self.wallx,self.wally,'b')             
        plt.plot(self.pointx,self.pointy,'r*')        
        plt.show()
        return 0
        
        
def main():
    
    des = Nozzle_Design()
    des.inishape()
    des.caculation()
    des.wall_caculation()
    
    des.display()


    
    return 0
    

if __name__=="__main__":
    main()