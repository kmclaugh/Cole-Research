def orbits_csv_to_comet_animation(orbit_filename,force_x_filename,force_y_filename):
    
    
    print("scanning position values")
    return_list = scan_position_values(orbit_filename)
    position_x_values = return_list[0]
    position_y_values = return_list[1]
    time_values = return_list[2]
    
       

    print("finding requested time values")
    counter = find_time_value(desired_time=0,time_list=time_values)

    print("scanning vector values")
    ## x-values
    vector_values = scan_vector_values(force_x_filename,counter)
    columbic_x_values = vector_values[0]
    radiation_x_values = vector_values[1]
    applied_wave_x_values = vector_values[2]
    velocity_x_values = vector_values[3]

    ## y-values
    vector_values = scan_vector_values(force_y_filename,counter)
    columbic_y_values = vector_values[0]
    radiation_y_values = vector_values[1]
    applied_wave_y_values = vector_values[2]
    velocity_y_values = vector_values[3]


    position_x_values = position_x_values[counter:]
    position_y_values = position_y_values[counter:]
    time_values = time_values[counter:]

    print("scaling vectors")
    ## columbic
    columbic_scale_factor = find_scale_factor(columbic_y_values,position_y_values)
    columbic_scaled_list_y = apply_scale_factor(columbic_y_values,columbic_scale_factor)
    columbic_scaled_list_x = apply_scale_factor(columbic_x_values,columbic_scale_factor)
    ## radiation
    radiation_scale_factor = find_scale_factor(radiation_y_values,position_y_values)
    radiation_scaled_list_y = apply_scale_factor(radiation_y_values,radiation_scale_factor)
    radiation_scaled_list_x = apply_scale_factor(radiation_x_values,radiation_scale_factor)
    ## applied wave
    applied_wave_scale_factor = find_scale_factor(applied_wave_y_values,position_y_values)
    applied_wave_scaled_list_y = apply_scale_factor(applied_wave_y_values,applied_wave_scale_factor)
    applied_wave_scaled_list_x = apply_scale_factor(applied_wave_x_values,applied_wave_scale_factor)
    ## velocity
    velocity_scale_factor = find_scale_factor(velocity_y_values,position_y_values)
    velocity_scaled_list_y = apply_scale_factor(velocity_y_values,velocity_scale_factor)
    velocity_scaled_list_x = apply_scale_factor(velocity_x_values,velocity_scale_factor)

    ## axis factors
    minx = min(position_x_values)
    miny = min(position_y_values)
    maxx = max(position_x_values)
    maxy = max(position_y_values)
    deltay = (maxy - miny)/12
    deltax = (maxx - minx)/12

    print("deleting old files")
    for the_file in os.listdir(animation_directory):
        file_path = os.path.join(animation_directory, the_file)
        try:
            os.unlink(file_path)
        except Exception as e:
            pass
    
    print("making animations files")
    counter = 1
    while counter <= len(time_values):
        tracker = log10(counter)
        if tracker % 1 == 0:
            print("creating file: {}".format(counter))   
        ## Values to be plotted for this frame
        plotx = []
        ploty = []
        if counter < 4:
            ##get values for this frame
            plotx += position_x_values[:counter]
            ploty += position_y_values[:counter]
            t = time_values[counter-1]
            ##create the time textbox
            txt1 = "Time: {} s".format(t)
            txt1 = txt1[:8]+txt1[17:]
            text(minx, maxy, txt1, fontsize=14,verticalalignment='top')
            ##Define the axis
            axis([minx-deltax,maxx+deltax, miny-deltay,maxy+deltay])
            ## Calculate Vector info
            ##columbic vector
            arrow(plotx[counter-1],ploty[counter-1],columbic_scaled_list_x[counter-1],columbic_scaled_list_y[counter-1],fc="w",ec="g",head_width=.1*deltay,head_length=.1*deltay, width=.005*deltay)
            ##radiation vector
            arrow(plotx[counter-1],ploty[counter-1],radiation_scaled_list_x[counter-1],radiation_scaled_list_y[counter-1],fc="w",ec="k",head_width=.1*deltay,head_length=.1*deltay, width=.005*deltay)
            ##applied wave vector
            arrow(plotx[counter-1],ploty[counter-1],applied_wave_scaled_list_x[counter-1],applied_wave_scaled_list_y[counter-1],fc="w",ec="m",head_width=.1*deltay,head_length=.1*deltay, width=.005*deltay)
            ##velocity vector
            arrow(plotx[counter-1],ploty[counter-1],velocity_scaled_list_x[counter-1],velocity_scaled_list_y[counter-1],fc="w",ec="r",head_width=.1*deltay,head_length=.1*deltay, width=.005*deltay)
     
            ##Plot tail
            plot(plotx,ploty,marker=None, linestyle="-", color="b")
            ##Plot particle
            plot(plotx[counter-1],ploty[counter-1],marker="o", linestyle="-", color="b")    
            ##Save figure
            filename = animation_directory+("/orbit_animation_%010d.png" % counter)
            savefig(filename, bbox_inches=0)
            clf()
        else:
            count_down = counter-4
            for i in range(count_down,counter):
                plotx.append(position_x_values[i])
                ploty.append(position_y_values[i])
           
            
            
            ##Define the axis
            axis([minx-deltax,maxx+deltax, miny-deltay,maxy+deltay])

            ## Calculate Vector info
            ##columbic vector
            arrow(plotx[-1],ploty[-1],columbic_scaled_list_x[i],columbic_scaled_list_y[i],fc="w",ec="g",head_width=.1*deltay,head_length=.1*deltay, width=.005*deltay)
            ##radiation vector
            arrow(plotx[-1],ploty[-1],radiation_scaled_list_x[i],radiation_scaled_list_y[i],fc="w",ec="k",head_width=.1*deltay,head_length=.1*deltay, width=.005*deltay)
            ##applied wave vector
            arrow(plotx[-1],ploty[-1],applied_wave_scaled_list_x[i],applied_wave_scaled_list_y[i],fc="w",ec="m",head_width=.1*deltay,head_length=.1*deltay, width=.005*deltay)
            ##velocity vector
            arrow(plotx[-1],ploty[-1],velocity_scaled_list_x[i],velocity_scaled_list_y[i],fc="w",ec="r",head_width=.1*deltay,head_length=.1*deltay, width=.005*deltay)

            ##create the time textbox
            t = time_values[i]
            txt = "Time: {} s".format(t)
            txt = txt[:8]+txt[17:]
            text(minx, maxy, txt, fontsize=14,verticalalignment='top')
            ##plot
            plot(plotx,ploty,marker=None, linestyle="-", color="b")
            plot(plotx[-1],ploty[-1],marker="o", linestyle="-", color="b")
            filename = animation_directory+("/orbit_animation_%010d.png" % counter)
            savefig(filename, bbox_inches=0)
            clf()
        counter += 1
        if counter > 300:
            break
    print("stitching animation")
    bash_command = "~/ffmpeg/ffmpeg -r 5 -qscale 0 -i orbit_animation_%10d.png ~/Desktop/orbit_animation{}.mp4".format(datetime.datetime.now())
    os.system(bash_command)
    textsender(8173126800,"Sir, simulation complete is complete. There were {} files created".format(counter))
