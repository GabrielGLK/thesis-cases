#! /usr/bin/python
# -*- coding: utf-8 -*-

# this file use to draw the images to explain
# how to get the volume of the bubbles

from __future__ import unicode_literals
import os
import sys
import numpy as np
import math
import matplotlib 
import matplotlib.pyplot as plt
import pylab
import argparse
import string
import math
import re
import operator
from scipy import ndimage


from PIL import Image
from PIL import ImageEnhance
from matplotlib.backends.backend_pdf import PdfPages


#from DataP import DataProcess
#import TaylorCorrelation as TC
#import lookup as LOOKUP

matplotlib.rcParams['text.usetex'] = True
matplotlib.rcParams['text.latex.unicode'] = True
matplotlib.rcParams['font.family'] = 'serif'
matplotlib.rcParams['font.size'] = 10
matplotlib.rcParams['xtick.direction'] = 'in'
matplotlib.rcParams['ytick.direction'] = 'in'
matplotlib.rcParams['figure.dpi'] = 50
matplotlib.rcParams['figure.dpi'] = 50

CWP      = os.path.abspath("./")
dir_list = [60,70,80,90,95]
case = [1,1.12,1.24,1.33,1.45,1.6,1.72,2,2.2,2.5,2.8,3]
a = 0
i = 6
while i < 12:  
    try:
        PATH_REFERENCE = os.path.abspath('/media/lguo/Elements1/calculation/Taylor_bubble/2D/singularity/sudden/expansion/ratios/wide_range/wide_range/expansion_wide/{}/{}'.format(str(dir_list[a]),str(case[i])))
        PATH_REFERENCE_1 = os.path.abspath('/media/lguo/Elements/calculation/test-cases/no-phase-change/taylor-bubble/2D/sudden_expansion_same/60/4/taylor_2D/velocity')
        PATH_INTERFACE = PATH_REFERENCE + '/interface'
        PATH_VELOCITY =  PATH_REFERENCE + '/velocity'
        PATH_OUT = PATH_REFERENCE + '/log_1'
        #PATH_COM = os.path.abspath("/media/lguo/Elements1/writing/Article/elsarticle-taylor/figure/2D/velocity_vector/straight")
        PATH_COM = os.path.abspath("/media/lguo/Elements1/calculation/Taylor_bubble/2D/singularity/sudden/expansion/ratios/wide_range/wide_range/expansion_wide/figure")
        PATH_FIG = PATH_COM + '/{}/{}-{}'.format(str(dir_list[a]),str(dir_list[a]),str(case[i]))
        PATH_FIG_SAVE = os.path.abspath('/media/lguo/Elements/phd-manuscript/thesis-manuscript/Taylor-bubble/latex/elsarticle-taylor/figure/mesh/mesh')
        if not os.path.exists(PATH_COM):
            os.makedirs(PATH_COM)
        if not os.path.exists(PATH_FIG):
            os.makedirs(PATH_FIG)

        def get_sur(filename):
            fn = os.path.join(PATH_REFERENCE_1 , filename)
            sur = []
            with open(fn, 'r') as file:
                seg = []
                for line in file.readlines():
                    line = line.rstrip('\n')
                    l = line.split(' ')
                    if len(l) == 2:
                        #print float(l[0]), float(l[1])
                        seg.append([float(l[0]), float(l[1])])
                    else:
                        segn = seg[:]
                        sur.append(segn)
                        seg = []
            return sur

        def get_out(filename):
            fn = os.path.join(PATH_OUT, filename)
            sur = []
            with open(fn, 'r') as file:
                seg = []
                for line in file.readlines():
                    line = line.rstrip('\n')
                    l = line.split(' ')
                    if len(l) == 3:
                        #print float(l[0]), float(l[1])
                        seg.append([float(l[0]), float(l[1]),float(l[2])])
                    else:
                        segn = seg[:]
                        sur.append(segn)
                        seg = []
            return sur

        def get_sur_1(time):
            fn = os.path.join(PATH_INTERFACE, 'interface-%g' % time)
            sur = []
            with open(fn, 'r') as file:
                seg = []
                for line in file.readlines():
                    line = line.rstrip('\n')
                    l = line.split(' ')
                    if len(l) == 2:
                        #print float(l[0]), float(l[1])
                        seg.append([float(l[0]), float(l[1])])
                    else:
                        segn = seg[:]
                        sur.append(segn)
                        seg = []
            return sur

        def get_sur_reference(filename):
            fn = os.path.join(PATH_REFERENCE, filename)
            arrs = []
            strr = ''
            with open(fn, 'r') as file:
                for line in file:
                    if '1e+30' in line:
                        line=line.replace('1e+30','0')
                    strr += line
            with open(fn, 'w') as file:
                file.write(strr)
            with open(fn, 'r') as file:
                for line in file.readlines():
                    line = line.strip().strip('\n')
                    l = line.split(' ')
                    #print float(l[0]), float(l[1]),float(l[2]),float(l[3])
                    arrs.append([float(l[0])-8, float(l[1]),float(l[2]),float(l[3])])
            return arrs

        def get_sur_2(time):
            fn = os.path.join(PATH_VELOCITY, 'vprof-%.1f' % time)
            arrs = []
            strr = ''
            with open(fn, 'r') as file:
                for line in file:
                    if '1e+30' in line:
                        line=line.replace('1e+30','0')
                    strr += line
            with open(fn, 'w') as file:
                file.write(strr)
            with open(fn, 'r') as file:
                for line in file.readlines():
                    line = line.strip().strip('\n')
                    l = line.split(' ')
                    #print float(l[0]), float(l[1]),float(l[2]),float(l[3])
                    arrs.append([float(l[0])-8, float(l[1]),float(l[2]),float(l[3])])
            return arrs

        def get_bubble_head(sur):
            maxx = sur[0][0][0]
            for seg in sur:
                x1 = seg[0][0]
                y1 = seg[0][1]
                x2 = seg[1][0]
                y2 = seg[1][1]
                if x1 > maxx:
                    maxx = x1
                if x2 > maxx:
                    maxx = x2
            return maxx

        def get_bubble_tail(sur):
            minn = sur[0][0][0]
            for seg in sur:
                x1 = seg[0][0]
                y1 = seg[0][1]
                x2 = seg[1][0]
                y2 = seg[1][1]
                if x1 < minn:
                    minn = x1
                if x2 < minn:
                    minn = x2
            return minn


        def plot_sur(plt, sur, sig, color):
            
            for seg in sur:
                x1 = seg[0][0]
                x2 = seg[1][0]
                y1 = seg[0][1]
                y2 = seg[1][1]
                plt.plot([y1, y2], [x1 - sig, x2 - sig], color,lw = 3)
                res, = plt.plot([-y1, -y2], [x1- sig, x2 - sig], color, lw = 3)
            return res
        '''
        def context2array(sur):
            return np.matrix([map(float, re.split('\s+', ln.strip()))
                for ln in sur.splitlines() if ln.strip()])'''


        def plot_contour(x_dim, y_dim, x_steps, y_steps, scalar_field, v_min, v_max, levels=None):
            from matplotlib import cm
            x, y = np.mgrid[-x_dim/2:x_dim/2:x_steps*1j, -y_dim/2:y_dim/2:y_steps*1j]
            cs = plt.contourf(x, y, scalar_field, zorder=1, cmap=plt.cm.jet, extent=[-x_dim/2.0, x_dim/2.0, -y_dim/2.0, y_dim/2.0], vmin=v_min, vmax=v_max, levels=levels)
            plt.colorbar(cs)
            return cs.levels


        def plot_velocity(plt, sur):
        
            from matplotlib.colors import LogNorm

            soa =np.array(sur) 
            X,Y,U,V = zip(*soa)
            M  = np.hypot(U, V)
            ax = plt.gca()
            YN = np.array([-val for val in Y])
            ax.quiver(Y,X,V,U,M, units = 'xy',angles='xy',scale_units='xy', scale=20,cmap=plt.cm.jet,clim=(0,2.5))
            qq=ax.quiver(YN,X,V,U,M,units = 'xy', angles='xy',scale_units='xy',scale=20,cmap=plt.cm.jet,clim=(0,2.5))
            #cbar = plt.colorbar(qq,shrink=0.9,aspect=25)
            #cbar_ticks = np.linspace(0., 2., num=11, endpoint=True) 
            #cbar.set_ticks(cbar_ticks) 
            #cbar.ax.set_ylabel('$U^*$')


        def out_data(filename):
            fn = os.path.join(PATH_OUT , filename)
            seg = []
            with open(fn, 'r') as file:
                for line in file.readlines():
                    line = line.rstrip('\n')
                    l = line.split(' ')
                    if len(l) == 26:
                        #print float(l[0]), float(l[1])
                        seg.append([float(l[0]), float(l[5]), float(l[6])])
                    else:
                        seg = []
            return seg
            

        def plot_legend(plt,time, xmin, ymin, xmax, ymax):
            #import case as case
            str_undim = 'log(Mo) = -3.65\nEo =41.86\nNf = 134.93\n$\\rho_r=934$\n$\mu_r=2906$'
            plt.axes().text(xmin - 2.8, (ymin+0.5), 
                str_undim,
                bbox={'facecolor':'white'}, 
                fontname='monospace',
                )
            str_time = r'$t^*$ = %.2f' % (time-14.1) 
            plt.axes().text(xmin - 2.8, (ymin-0.1), 
                str_time,bbox={'facecolor':'white'}, 
                fontname='monospace',
                color='r'
                )

        def plot_streamline(plt,sur,a,time,ex,y_lim,y_max,head,tail):
            from scipy.interpolate  import griddata
            ax = plt.gca()
            soa =np.array(sur) 
            x,y,u,v = zip(*soa)
            nx, ny =1000, 1000
            yn = np.array([-val for val in y])
            pts = np.vstack((x, y)).T
            ptsi = np.vstack((x, yn)).T
            vals = np.vstack((u, v)).T
            # lower-parts
            xi = np.linspace(-4,8,1000)
            yi = np.linspace(0,ex[i]*0.5+0.08 ,1000)
            yni = np.linspace(-ex[i]*0.5-0.08,0,1000)
            ipts = np.vstack(a.ravel() for a in np.meshgrid(yi, xi)[::-1]).T
            iptsi= np.vstack(a.ravel() for a in np.meshgrid(yni, xi)[::-1]).T
            ivals = griddata(pts, vals, ipts, method='cubic')
            ivalsi = griddata(ptsi, vals, iptsi, method='cubic')
            ui, vi = ivals.T
            ui_i,vi_i = ivalsi.T
            ui.shape = vi.shape = (ny, nx)
            ui_i.shape = vi_i.shape = (ny, nx)
            time_1 = time*100 + 1
            b = int(time_1)
            # upper-parts

            c = ex[i]*0.5

            ### define the position and range of the streamlines
            stream_points_1 = np.array(zip(np.arange(0.5,0,-0.1), np.arange(head+0.1-8,tail-1-8,-.1))) # x>0 for lower tube
            stream_points_2 = np.array(zip(np.arange(c,0.55,-0.1), np.arange(0,4,.1)))  # x>0 for upper tube
            stream_points_3 = np.array(zip(np.arange(-0.5,0,0.1), np.arange(head+0.1-8,tail-1-8,-.1))) # x<0 for lower tube
            stream_points_4 = np.array(zip(np.arange(-c,-0.55,0.1), np.arange(0,4,.1))) # x<0 for upper tube

            # vortex near bubble tail
            stream_points_5 = np.array(zip(np.arange(-0.05,0,0.05), np.arange(tail-8.1,tail-8.15,-.05)))
            stream_points_6 = np.array(zip(np.arange(0.05,0,-0.05), np.arange(tail-8.1,tail-8.15,-.05)))
            # vortex inside bubble 
            stream_points_7 = np.array(zip(np.arange(0.2,0,-0.1), np.arange(tail-7.5,head-9.5,.1)))
            stream_points_8 = np.array(zip(np.arange(-0.2,0,0.1), np.arange(tail-7.5,head-9.5,.1)))



            #### draw streamlines
            ax.streamplot(yi, xi, vi-a[b][2], ui-a[b][1], color = 'k',start_points=stream_points_1, arrowsize = 0.3,linewidth = 0.5,density =20)
            ax.streamplot(yi, xi, vi-a[b][2], ui-a[b][1], color = 'k',start_points=stream_points_2, arrowsize = 0.3,linewidth = 0.5,density =20)

            ax.streamplot(yni, xi, a[b][2] - vi_i, ui_i-a[b][1], color = 'k',start_points=stream_points_3, arrowsize = 0.3,linewidth = 0.5,density =20,integration_direction='both')
            ax.streamplot(yni, xi, a[b][2] - vi_i, ui_i-a[b][1], color = 'k',start_points=stream_points_4, arrowsize = 0.3,linewidth = 0.5,density =20,integration_direction='forward')

            ax.streamplot(yi, xi, vi-a[b][2], ui-a[b][1], color = 'k',start_points=stream_points_6, arrowsize = 0.3,linewidth = 0.5,density=20,integration_direction='both')
            ax.streamplot(yni, xi, a[b][2] - vi_i, ui_i-a[b][1], color = 'k',start_points=stream_points_5, arrowsize = 0.3,linewidth = 0.5,density =20,integration_direction='both')

            ax.streamplot(yi, xi, vi-a[b][2], ui-a[b][1], color = 'k',start_points=stream_points_7, arrowsize = 0.3,linewidth = 0.5,density =20,integration_direction='both')
            ax.streamplot(yni, xi, a[b][2] - vi_i, ui_i-a[b][1], color = 'k',start_points=stream_points_8, arrowsize = 0.3,linewidth = 0.5,density =20,integration_direction='both')
            
            
        def simulation(a):
       
            time_start = {'60':11.6,'70':12.7,'80':14.2,'90':12,'95':12}
            #time_end = {'60':sixty[str(case[i])],'70':seventy[str(case[i])],'80':eighty[str(case[i])],'90':ninty[str(case[i])],'95':nin_five[str(case[i])]}
            time_end = {'60':30,'70':30,'80':30,'90':25,'95':25}

            for time in np.arange(time_start[str(dir_list[a])],time_end[str(dir_list[a])],1):

                plt.figure(figsize=(4, 4))

                """
                Set labels
                """
                #plt.xlabel(r'$r/D$')
                #plt.ylabel(r'$(z_h - z)/D$')
                #plt.axis('off')
                """
                Set range
                """
                #plt.xscale('log')
                #plt.yscale('log')
                # plt.ylim([0,0.5])

                sur1 = get_sur_1(time)
                sur2 = get_sur_2(time)
                #plot_velocity(plt,sur2)
                a = out_data('out')

                ex= [1,1.12,1.24,1.33,1.45,1.6,1.72,2,2.2,2.5,2.8,3] # expansion ratio
                
                # change interface name to float
                '''filename_list = os.listdir(CWP)  
                for filename in filename_list:   
                    strlen = len(reply.encode(filename))
                    str_reverse = filename[::-1]  
                    start = str.find(filename,"-") + 1 
                    end = str.find(str_reverse,"(") 
                    end = strlen - end - 1
                    apxstart = str.find(str_reverse,".") 
                    apxstart = - apxstart - 1 
                    final_filename = filename[start:end] + filename[apxstart:]
                    print('interface-'+filename[start:end])   
                    print('interface=-1-:' + final_filename)'''
        


                
                xh1 = get_bubble_head(sur1)
                xh11 = get_bubble_tail(sur1)
                with open('data', 'wa') as f:  
                    f.write("head-%.1f :%f\t" %(time, xh1))
                    f.write("tail-%.1f :%f\n" %(time, xh11))
                le1 = plot_sur(plt, sur1, 8, "r")
                
                
                ex_ratio = ex[i]  
                #plt.title('$t = %.2f s$' %(0.04*(time- 14.2) ))
                if ex_ratio in ex[:]:
                    expansion = ex_ratio + 0.08
                    plt.xlim(-0.5*expansion,0.5*expansion)
                else:
                    expansion = 1 + 0.08
                    plt.xlim(-0.5*expansion,0.5*expansion)
    
                y_min = -3.5 + (time - 11.6)/3
                y_max = xh1-8+0.5
                #plot_streamline(plt,sur2,a,time,ex,y_min,y_max,xh1,xh11)

                plt.ylim(y_min,y_max)
                if y_min<0:
                    y1 = math.fabs(y_min)/(y_max - y_min)
                else:
                    y1 = 0
                x1 = math.fabs(0.5*expansion+0.5)/expansion
                x2 = math.fabs(0.5+expansion)/expansion - 0.08
                x3 = 0
                x4 = (0.5*expansion-0.5 )/expansion
                
                plt.axvline(x = 0.5,ymin = 0, ymax = y1, color='b')
                plt.axvline(x = -0.5,ymin = 0, ymax = y1, color='b')
                plt.axvline(x = 0.5*ex_ratio,ymin = y1, ymax = 1, color='b')
                plt.axvline(x = -0.5*ex_ratio,ymin = y1, ymax = 1, color='b')
                plt.axhline(y = 0,xmin=0.08*0.5/expansion,xmax=(expansion*0.5-0.5)/expansion,color='b')
                plt.axhline(y = 0,xmin=(0.5+0.5*expansion)/expansion,xmax=1-0.08*0.5/expansion,color='b')
                #plot_legend(plt, time, -0.5, y_min, 0.5, y_max)
                plt.axhspan(xmin = x1,xmax = x2,ymin = y_min,ymax = 0,facecolor = 'w',alpha = 1)
                plt.axhspan(xmin = x3,xmax = x4,ymin = y_min,ymax = 0,facecolor = 'w',alpha = 1)


            #plt.legend([le1,le2,le3], ["$h = 1/64$","$h = 1/128$","$h = 1/256$"], loc="upper right", bbox_to_anchor=(-0.3,1))

            #plt.axhline(xh11-xh1,xmin=0.7,xmax=1.2,color='b')
            #plt.annotate('$z_{head}$',xy=(0,0),xytext=(-0.5,0),weight='bold',color='r',arrowprops=dict(arrowstyle='->',connectionstyle='arc3',color='k'))
            #plt.annotate('$z_{tail}$',xy=(0.35,xh33-xh3),xytext=(-0.2,xh33-xh3-0.05),weight='bold',color='r',arrowprops=dict(arrowstyle='->',connectionstyle='arc3',color='k'))

            # plt.grid(True)
                version = matplotlib.__version__

                plt.axes().set_aspect('equal')
                plt.tight_layout()

                name =  "basilisk_" + "%g"% (4*(time- 11.6))
                #fn = os.path.join(PATH_FIG, (name + ".png"))
                fn = os.path.join(PATH_FIG, (name + ".pdf"))
                pdf = PdfPages(fn.format(version))
                #plt.savefig(file_path + '.png', dpi=Vc.dpi)
                #fig = plt.figure()
                #fig.savefig(fn, dpi = fig.dpi)
                plt.savefig(fn, format='pdf', bbox_inches='tight')
                plt.close() 
            
            # plt.show()


        def main():
            simulation(a)


        if __name__ == '__main__':
            print("-----start-----")
            main()  
            print("-----finished!!!-----")
    except ValueError:
        pass
    i += 1