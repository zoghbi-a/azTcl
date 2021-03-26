puts "\n@@ loading az tools @@\n"

# Use steppar to scan the residuals of a fit to search for lines
proc az_scan_en_norm { args } {
    # Use steppar to scan the residuals of a fit
    # to search for lines

    # print help if requested or if no arguments are given #
    if {![llength $args] || [lindex $args 0] == "?" || [lindex $args 0] == "-h"} {
        set help "\n\
    Use steppar to scan the residuals of a fit\n\

    Parameters:\n\
        xcm_file:  name of of file to load data and model.\n\
        inew_comp: where to add the searching zgauss line\n\
        en_pars:   energy parameters of steppar {ipar start stop nstep}\n\
        norm_pars: norm parameters of steppar {ipar start stop nstep}\n\
        mode:      default:0; if 1: do opposite sign too
        redshift:  redshift for the scanning gaussian
        "   
        puts $help
        return
    }


    # check number of arguments #
    if {[llength $args] != 6} {
        puts "\nI need 6 arguments. See -h for help. Quitting ...\n"
        return
    }

    # unpack the arguments #
    set xcm_file  [lindex $args 0]
    set inew_comp [lindex $args 1]
    set en_pars   [lindex $args 2]
    set norm_pars [lindex $args 3]
    set mode      [lindex $args 4]
    set z         [lindex $args 5]


    # log? #
    set p0 [lindex $en_pars 0]
    set en_log "nolog"
    if {$p0 == "log"} {set en_log "log"}
    if {$p0 == "log" || $p0 == "nolog"} { set en_pars [lrange $en_pars 1 4] }

    set p0 [lindex $norm_pars 0]
    set norm_log "nolog"
    if {$p0 == "log"} {set norm_log "log"}
    if {$p0 == "log" || $p0 == "nolog"} { set norm_pars [lrange $norm_pars 1 4] }



    # output file in xcm_file[-xcm].scan
    set out_file [format "%s.scan" [string range $xcm_file 0 \
            [expr [string length $xcm_file]-5]]]
    

    # mode #
    if {$mode != 0 && $mode != 1} {
        puts "mode need to be 0 or 1; Quitting ..."
        return
    }

    # save query and chatter values, so reset them at the end #
    set query [tcloutr query]
    set chatlevel [scan [tcloutr chatter] "%d"]
    query yes
    #chatter 0
    ## ---------------------- ##

    

    ## -- load xcm file -- ##
    @$xcm_file
    fit 500
    set stat0 [tcloutr stat]
    scan [tcloutr dof] "%d %d" dof npha
    ## ------------------- ##

    
    ## -- Starting steppar -- ##
    puts "@@@ Starting steppar in energy and norm @@@"
    addcomp $inew_comp zgauss & 6.4 & 0 -0.01& $z& 0.0 0.01 -1e6 -1e6 1e6 1e6 & /*
    puts "steppar best $en_log $en_pars $norm_log $norm_pars"
    steppar best $en_log $en_pars $norm_log $norm_pars
    set norm [tcloutr steppar [lindex $norm_pars 0]]
    set en   [tcloutr steppar [lindex $en_pars 0]]
    set stat [tcloutr steppar stat]

    if {$mode == 1} {
        puts "@@ doing opposite sign @@"
        add $inew_comp const & -1 -1 -1 -1& /*
        lset en_pars 0 [expr [lindex $en_pars 0]+1]
        lset norm_pars 0 [expr [lindex $norm_pars 0]+1]
        puts "steppar best $en_log $en_pars $norm_log $norm_pars"
        steppar best $en_log $en_pars $norm_log $norm_pars
        set norm_neg [tcloutr steppar [lindex $norm_pars 0]]
        set en_neg   [tcloutr steppar [lindex $en_pars 0]]
        set stat_neg [tcloutr steppar stat]
        delcomp $inew_comp
    }
    delcomp $inew_comp
    ## ---------------------- ##


    # write the output to a file
    # values are not sorted. it is easier to sort them outside tcl
    puts "@@@ Writing the results to $out_file @@@"
    set nval [llength $en]
    set fp [open $out_file w]
    puts $fp "# chi0: $stat0 ; dof: $dof ; npha: $npha"
    
    for {set ic 0} {$ic < $nval} {incr ic} {
        puts -nonewline $fp [format "%.5g %.5g %.5g\n" \
            [lindex $en $ic] [lindex $norm $ic] [lindex $stat $ic]]
    }

    # write the opposite sign if requested #
    if {$mode == 1} {
        for {set ic 0} {$ic < $nval} {incr ic} {
            puts -nonewline $fp [format "%.5g %.5g %.5g\n" \
            [lindex $en_neg $ic] -[lindex $norm_neg $ic] [lindex $stat_neg $ic]]
        } 
    }

    close $fp
    query $query
    chatter $chatlevel
}


proc az_sim_dchi2 { args } {
    # Simulate dchi2 improvement due to a narrow gaussian line
    # for some fit to the data

    # print help if requested or if no arguments are given #
    if {![llength $args] || [lindex $args 0] == "?" || [lindex $args 0] == "-h"} {
        set help "\n\
    Simulate dchi2 improvements due to a narrow gaussian line\n\

    Parameters:\n\
        niter: number of simulations to run.\n\
        xcm_file:  name of of file to load data and model.\n\
        inew_comp: where to add the new gauss line. e.g. 1 --> addcomp 1 gauss\n\
        outfile:   energy parameters of steppar {ipar start stop nstep}. 0(default): xcm->sim\n\
        "   
        puts $help
        return
    }


    # check number of arguments #
    if {[llength $args] < 3} {
        puts "\nI need 6 arguments. See -h for help. Quitting ...\n"
        return
    }

    # unpack the arguments #
    set niter     [lindex $args 0]
    set xcm_file  [lindex $args 1]
    set inew_comp [lindex $args 2]
    if {[llength $args] == 4} {
        set outfile [lindex $args 3]
    } else {
        set outfile 0
    }
    

    # file names #
    set base [file rootname $xcm_file]
    if {$outfile == 0 } { set outfile $base.sim }
    set modfile ${base}_mod.xcm


    ## -- chatter & qury -- ##
    scan [tcloutr chatter] "%d" chatter
    chatter 0
    set query [tcloutr query]
    query yes


    ## -- read xcm file -- ##
    @$xcm_file
    fit 1000

    ## -- save model -- ##
    rm $modfile &> /dev/null
    save model $modfile
    

    ## -- save filenames -- ##
    set ndata [tcloutr datasets]
    for {set id 1} {$id <= $ndata} {incr id} {
        set sfile($id) [tcloutr filename $id]
        set bfile($id) [string trim [tcloutr backgrnd $id]]
        if {$bfile($id) eq ""} {set bfile($id) "none"}
        set rspfile($id) [string trim [tcloutr response $id]]
        if {$rspfile($id) eq ""} {set rspfile($id) "none"}
        set arffile($id) [string trim [tcloutr arf $id]]
        if {$arffile($id) eq ""} {set arffile($id) "none"}
        # save the noticed parameters
        set noticed($id) [tcloutr noticed $id]
    }
    ## -------------------- ##


    ## -- save parameters -- ##
    set npar [tcloutr modpar]
    for {set ip 1} {$ip <= $npar} {incr ip} {
        scan [tcloutr param $ip] "%f" pval($ip)
        set pislink($ip) [string range [tcloutr plink $ip] 0 0]
        set plink($ip) [string trimleft [tcloutr plink $ip] ?TF?]
    }
    ## --------------------- ##

    ## -- original stat -- ##
    set stat0 [tcloutr stat]
    ## ------------------- ##

    #--- write to output ---#
    set fileid [open $outfile w]
    puts $fileid [format "# chi0: %6.6g" $stat0]
    close $fileid
    #-----------------------#


    ## --------------------- ##
    ## -- Start main loop -- ##
    ## --------------------- ##
    for {set iter 1} { $iter<=$niter } {incr iter} {

        puts "@@@@@ iter $iter @@@@@"

        #-- get the random set of parameters --#
        set simp [tcloutr simpars]
        set outstr ""
        for {set ip 0} {$ip<$npar} {incr ip} {
            set ii [ expr $ip+1 ]
            if { $pislink($ii) == "T" } {
                append outstr "& $plink($ii)"
            } else {
                append outstr "& $pval($ii)"
            }
        }
        newpar 1-$npar $outstr
        #--------------------------------------#

        #--- fake the spectra ---#
        set fakep ""
        for {set i 1} {$i<=$ndata} {incr i} {append fakep "&&&&&"}
        fakeit nowrite $fakep
        #------------------------#

        
        #--- ignore the appropriate channels ---#
        model none
        for {set iset 1} {$iset <= $ndata} {incr iset} {
            ignore $iset:**
            notice $iset:$noticed($iset)
        }
        #----------------------------------------#

        #--- load original model and fit ---#
        @$modfile
        fit; fit; fit
        set stat [tcloutr stat]
        #------------------------------------#

        #-- max residuals: +ve and -ve --#
        set resid [ tcloutr peakrsid]
        # +ve #
        addcomp $inew_comp gauss & [lindex $resid 0] & 0 -1&\
                 [lindex $resid 1] 0.01 -1e6 -1e6 1e6 1e6
        fit; fit ; fit
        set stat1 [tcloutr stat]
        delcomp $inew_comp
        # -ve #
        addcomp $inew_comp gauss & [lindex $resid 2] & 0 -1&\
                 [lindex $resid 3] 0.01 -1e6 -1e6 1e6 1e6
        fit; fit ; fit
        set stat2 [tcloutr stat]
        delcomp $inew_comp
        # ---- #
        if { $stat2 < $stat1 } { set stat1 $stat2 }
        set dstat [expr $stat - $stat1 ]
        #--------------------------------#

        #--- write to output ---#
        set fileid [open $outfile a]
        puts $fileid [format "%d %6.6g %6.6g %6.6g" $iter $stat $stat1 $dstat]
        close $fileid
        #-----------------------#

        
        #---------------#
        #---- RESET ----#
        #---------------#
        model none
         for {set iset 1} {$iset <= $ndata} {incr iset} {
            data [expr $ndata+$iset] $sfile($iset)
            response [expr $ndata+$iset] $rspfile($iset)
            if { $bfile($iset) ne "none" } {
                back [expr $ndata+$iset] $bfile($iset)
            }
            if { $arffile($iset) ne "none" } {
                arf [expr $ndata+$iset] $arffile($iset)
            }
        }
        for {set iset 1} {$iset <= $ndata} {incr iset} {
            data none /
        }

        # load original model #
        @$modfile
        fit
        #---------------#
    }
    query $query
    chatter $chatter

}


# get list of free params
proc az_free_params {args} {
    # get list of free params
    # print help if requested #
    if {[lindex $args 0] == "?" || [lindex $args 0] == "-h"} {
        set help "\n\
    get list of free params\n\
        "   
        puts $help
        return
    }

    set npar [tcloutr modpar]
    set varpar {}

    for {set ip 1} {$ip<=$npar} {incr ip} {
        if {[info exists frozen]} { unset frozen }
        scan [tcloutr param $ip] "%f %f" dum frozen
        scan [tcloutr plink $ip] "%s" linked
        if { ! [info exists frozen] } {continue}
        if { $linked == "T" || $frozen < 0 } {continue}
        lappend varpar $ip
    }
    return $varpar
}


# find parameter in multiple data groups 
proc az_find_ipar { args } {

    # print help if requested or if no arguments are given #
    if {![llength $args] || [lindex $args 0] == "?" || [lindex $args 0] == "-h"} {
        set help "\n\
    return the parameter number from multiple data groups \n\

    Parameters:\n\
        ip:  parameter number from the first group\n\
        "   
        puts $help
        return
    }

    # check number of arguments #
    if {[llength $args] > 1} {
        puts "\nI need at most 1 argument. See -h for help. Quitting ...\n"
        return
    }

    # unpack the arguments #
    set ip [lindex $args 0]

    # how many parameters per dataset do we have?
    set ngrp [tcloutr datagrp]
    set npar [expr [tcloutr modpar]/$ngrp ]

    # now untie the ipar's from datasets n #
    set text {}
    for {set i 0} {$i<$ngrp} {incr i} {
        lappend text [expr $i * $npar + $ip]
    }
    return $text
}



# Untie parameters of a spectrum
proc az_untie_spec { args } {
    # Untie parameters of a spectrum

    # print help if requested or if no arguments are given #
    if {![llength $args] || [lindex $args 0] == "?" || [lindex $args 0] == "-h"} {
        set help "\n\
    Untie parameters of a spectrum\n\

    Parameters:\n\
        n:  spectrum number (or a list) whose parameters are to be untied\n\
        ipar (optional): a list of parameters to untie. If not given, untie all\n\
        "   
        puts $help
        return
    }

    # check number of arguments #
    if {[llength $args] > 2} {
        puts "\nI need at most 2 arguments. See -h for help. Quitting ...\n"
        return
    }

    # unpack the arguments #
    set n [lindex $args 0]
    if {[llength $args] == 2} {
        set ipar [lindex $args 1]
    } else {
        set ipar {}
    }
    

    # how many parameters per dataset do we have?
    set npar [expr [tcloutr modpar]/[tcloutr datagrp] ]

    # make ipar if not given 
    if {[llength $ipar] == 0} {
        for {set i 1} {$i <= $npar} {incr i} {
            lappend ipar $i
        }
    }

    # now untie the ipar's from datasets n #
    set text ""
    for {set i 0} {$i<[llength $n]} {incr i} {
        foreach j $ipar {
            append text " [expr ([lindex $n $i]-1) * $npar + $j]"
        }
    }
    puts "untie $text"
    untie $text
}

# calculate error for a single parameter #
proc _az_error {delchi2 ipar {new_fit 1}} {
    # calculate error for a single parameter 
    # delchi2
    # ipar: parameter number
    # new_fit: 0: return even if a new fit is found.
    #          1: recalculate error after new fit


    while {1} {
        fit 100 1e-2
        catch { error max 1000 nonew [expr $delchi2*1.0] $ipar }
        scan [tcloutr error $ipar] {%f %f %s} minerr maxerr flags
        scan [tcloutr param $ipar] {%f %f %f %f %f %f} val dum hmin dum dum hmax

        if { ($new_fit == 0) || ( [string index $flags 0] == "F" ) } { 
            break 
        }
    }

    # convert from min/max to +/- errors #
    if { [string index $flags 3] == "T" } {set minerr $hmin}
    if { [string index $flags 4] == "T" } {set maxerr $hmax}
    if {$minerr==$val} {set minerr $hmin}
    if {$maxerr==$val} {set maxerr $hmax}
    set res [format "%.4e %.4e %.4e" \
        $val [expr $maxerr-$val] [expr $minerr-$val]]

    # if new_fit=0; append T|F if a new fit is found #
    if { $new_fit == 0 } {
        append res " [string index $flags 0]"
    }
    # val +err -err T|F (new fit found)
    return $res
}

# calculate error for a list of parameters #
proc _az_error_list {delchi2 {ipar_list 0} {reorder 0}} {
    # calculate error for a list of parameters
    # delchi2
    # ipar (optional): a list of parameter numbers. If not give
    #    use all free parameters (call az_free_params)
    # reorder: how to reorder parameters if a new fit is found
    #   0: icurr, {iprev_order}-icurr
    #   1: icurr, {iprev>icurr}, {iprev<icurr}
    
    set query [tcloutr query]
    scan [tcloutr chatt] "%d" chatt
    query yes
    chatter 0

    # list of parameters 
    if {[llength $ipar_list]==1 && $ipar_list <= 0} {
        set ipar_list [az_free_params]
    }
    set ipar_copy $ipar_list
    set npar [llength $ipar_list]
    set err_list {}

    # temporary xcm file #
    set tmpfile _az_[clock millisec].xcm
    fit 100 1e-2
    catch [file delete $tmpfile]
    save all $tmpfile


    # loop until we finish. If a new fit is found #
    # re-calculate errors starting from the current param
    set ic 0
    while {1} {
        incr ic
        set newfit 0
        for {set ip 0} {$ip<$npar} {incr ip} {

            set ipar [lindex $ipar_copy $ip]
            puts "\n-- errors for parameter $ipar :: [tcloutr pinf $ipar] ---"

            set err [_az_error $delchi2 $ipar 0]
            if {[lindex $err 3] == "T"} {
                set newfit 1
                puts "\n@@----------@@"
                puts "new best fit: [tcloutr stat]"
                # update ipar_copy, putting the current param first #
                set ipar_cp {}
                if {$ic < 10} { lappend ipar_cp $ipar }
                
                if {$reorder == 0} {
                    foreach iip $ipar_copy {
                        if {$iip != $ipar} {lappend ipar_cp $iip}
                    }
                } else {
                    foreach iip $ipar_copy {
                        if {$iip > $ipar} {lappend ipar_cp $iip}
                    }
                    foreach iip $ipar_copy {
                        if {$iip < $ipar} {lappend ipar_cp $iip}
                    }
                }
                if {$ic >= 10} { lappend ipar_cp $ipar; set ic 0}
                
                set ipar_copy $ipar_cp
                set err_list {}
                catch [file delete $tmpfile]
                save all $tmpfile
                puts "@ parameter list @"
                puts "$ipar_copy"
                puts "temporary xcm: $tmpfile"
                puts "@@----------@@\n"
                break
            }
            puts "   :: $err"
            lappend err_list [lrange $err 0 2]
        }
        if {$newfit == 0} {break}
    }

    # re-order err_list #
    set err_ord {}
    for {set i 0} {$i<$npar} {incr i} {
        for {set j 0} {$j<$npar} {incr j} {
            if {[lindex $ipar_list $i] == [lindex $ipar_copy $j]} {
                lappend err_ord [lindex $err_list $j]
            }
        }
    }
    query $query
    chatter $chatt
    return $err_ord
}

# calculate & print errors
proc az_calc_errors { args } {
    # calculate & print errors from a fit #

    # print help if requested or if no arguments are given #
    if {![llength $args] || [lindex $args 0] == "?" || [lindex $args 0] == "-h"} {
        set help "\n\
    Calculate & print errors from a fit\n\

    Parameters:\n\
        ipar: a list of parameter numbers\n\
        suff: a suffix for the output files: {suff}.xcm|log\n\
        delchi2 (optional): delta chi2 for the errors. Default 1.0 (for 1-sigma)\n\
        reorder (optional): how to reorder parameters when a new\
        fit is found.\
        0: icurr, {iprev_order}-icurr\
        1: icurr, {iprev>icurr}, {iprev<icurr}\n\
        "   
        puts $help
        return
    }

    # check number of arguments #
    if { [llength $args] < 2 || [llength $args] > 4 } {
        puts "\nWrong number of arguments. See -h for help. Quitting ...\n"
        return
    }

    # unpack the arguments #
    cpd none
    set ipar [lindex $args 0]
    set suff [lindex $args 1]
    set delchi2 1.0
    if {[llength $args] > 2} {
        set delchi2 [lindex $args 2]
    }
    set reorder 0
    if {[llength $args] > 3} {
        set reorder [lindex $args 3]
    }

    # get a error list #
    set err_list [_az_error_list $delchi2 $ipar $reorder]

    # fit statistics #
    set chi2 [tcloutr stat]
    scan [tcloutr dof] "%d %d" dof dum
    set rchi2 [expr 1.*$chi2/$dof]
    set nullp -1
    if {[tcloutr statmethod] == "Chi-Squared"} {
        set nullp [tcloutr nullhyp]
    }

    # save xcm file #
    rm $suff.xcm >& /dev/null
    save all $suff
    puts "-- saved $suff.xcm --"


    # write log file #
    set text [format "\n# %g %d %g %g\n" $chi2 $dof $rchi2 $nullp]
    for {set ip 0} {$ip<[llength $ipar]} {incr ip} {
        set name [ tcloutr pinfo [lindex $ipar $ip] ]
        set err [lindex $err_list $ip]
        set err_avg [expr ([lindex $err 1] - [lindex $err 2])/2]
        set err [linsert $err 1 [format "%.4e" $err_avg]]
        append text "$err \"$name\"\n"
    }
    set fp [open $suff.log w]
    puts $fp $text
    close $fp
    puts "-- saved $suff.log --"
}


# randomize parameters using tcloutr sim
proc az_rand_pars {args} {
    
    # print help if requested or if no arguments are given #
    if {[lindex $args 0] == "?" || [lindex $args 0] == "-h"} {
        set help "\n\
    Randomize parameters using either fit covariance or loaded mcmc chain\n\
    
    Parameters:\n\
        ret (optional): 0|1; return the newpar string without calling newpar. default 1\n\
        "   
        puts $help
        return
    }
    
    # process arguments 
    set ret 0
    if {[llength $args] > 0} {set ret [lindex $args 0]}


    # generate random parameters #
    set rand_pars [tcloutr sim]
    set free [az_free_params]
    set npar [tcloutr modpar]
    
    set ic 0
    set ifree [lindex $free $ic]
    set text "&"
    for {set i 1} {$i <= $npar} {incr i} {
        if {$i == $ifree} {
            append text [lindex $rand_pars [expr $i-1]]
            incr ic
            set ifree [lindex $free $ic]
        }
        append text "&"
    }
    if {$ret == 1} { return $text }
    newpar 1-$npar $text
}


# plot randomized model 
proc az_rand_plot {args} {
    # plot randomized model  #

    # print help if requested or if no arguments are given #
    if {![llength $args] || [lindex $args 0] == "?" || [lindex $args 0] == "-h"} {
        set help "\n\
    Plot randomized model\n\

    Parameters:\n\
        nsim: the number simulations\n\
        outfile: name of the output file\n\
        ptype (optional): plot type: mo, em, eem; default mo\n\
        igrp (optional): data group; default 1\n\
        rspcmd (optional): dummyrsp argument of use dummy response. e.g. \"2. 10. 100 log\"
        "   
        puts $help
        return
    }

    # check number of arguments #
    if { [llength $args] < 2 || [llength $args] > 5 } {
        puts "\nWrong number of arguments. See -h for help. Quitting ...\n"
        return
    }

    # unpack the arguments #
    set nsim    [lindex $args 0]
    set outfile [lindex $args 1]
    set ptype mo
    if {[llength $args] > 2} {set ptype [lindex $args 2]}
    set igrp 1
    if {[llength $args] > 3} {set igrp [lindex $args 3]}
    set rspcmd 0
    if {[llength $args] > 4} {set rspcmd [lindex $args 4]}


    scan [tcloutr chatt] "%d" chatt
    chatter 0

    # generate random parameters #
    set rand_pars {}
    for {set isim 1} {$isim<=$nsim} {incr isim} {
        lappend rand_pars [tcloutr sim]
    }

    # add a dummy response
    if {$rspcmd != 0} { 
        da none
        dummyrsp $rspcmd
    }


    # x-axis & free parameter indexes #
    set en   [tcloutr plot $ptype x $igrp]
    set ipar [az_free_params]


    # now do the simulations #
    set text ""
    for {set ie 0} {$ie < [llength $en]} {incr ie} {
        append text [format " %.6e" [lindex $en $ie]]
    }
    for {set isim 1} {$isim<=$nsim} {incr isim} {
        set sim [lindex $rand_pars [expr $isim-1]]
        set newp "1-[tcloutr modpar] &"
        for {set ip 1} {$ip<=[tcloutr modpar]} {incr ip} {
            if {[lsearch -exact $ipar $ip] >= 0} {
                append newp "[lindex $sim [expr $ip-1]] &"
            } else {
                append newp "&"
            }
        }
        newp $newp
        set plt [tcloutr plot $ptype y $igrp]
        append text "\n"
        for {set ie 0} {$ie < [llength $en]} {incr ie} {
            append text [format " %.6e" [lindex $plt $ie]]
        }
    }
    set fp [open $outfile.tmp w]
    puts $fp $text
    close $fp

    # use python to organize the output #
    set py "import numpy as np\n"
    append py "d = np.loadtxt('$outfile.tmp')\n"
    append py "en,mod = d\[0\], d\[1:\]\n"
    append py "mm, ms = mod.mean(0), mod.std(0)\n"
    append py "d = np.vstack((en,mm,ms,mod)).T\n"
    append py "txt = 'descriptor en mod,+- %s\\n'%"
    append py " ' '.join(\['s%d'%(i+1) for i in range(len(mod))\])\n"
    append py "txt += '\\n'.join(\[' '.join(\['%g'%s for s in x\]) for x in d\])\n"
    append py "with open('$outfile', 'w') as fp: fp.write(txt)\n"

    set fp [open $outfile.tmp.py w]
    puts $fp $py
    close $fp
    python $outfile.tmp.py
    file delete $outfile.tmp
    file delete $outfile.tmp.py

    chatter $chatt
}

# save current parameters to a string #
proc az_save_pars {} {
    # return a string that has the current parameters
    # so that the params can be loaded again by doing
    # newpar $res
    set n [tcloutr modpar]
    set txt " 1-$n "
    for {set ip 1} {$ip<=$n} {incr ip} {
        set link [tcloutr plink $ip]
        if { [string index $link 0] == "T" } {
            set lp [string trim $link "T = "]
            append txt "&=$lp"
        } else {
            append txt "&" [tcloutr param $ip]
        }
    }
    return $txt
}


# set parameters to best in chain #
proc az_best_chain_pars {} {
    # return a string  variable than when used as newpar $txt sets
    # the parameters to the best values in the currenlty loaded chain
    set b [tcloutr chain best]
    set ii [az_free]
    set ic 0
    
    set n [tcloutr modpar]
    set txt " 1-$n "
    for {set ip 1} {$ip<=$n} {incr ip} {
        set link [tcloutr plink $ip]
        if { [string index $link 0] == "T" } {
            set lp [string trim $link "T = "]
            append txt "&=$lp"
        } elseif { $ic < [llength $ii] && $ip == [lindex $ii $ic] } {
            append txt "& " [lindex $b $ic]
            incr ic
        } else {
            append txt "&" [tcloutr param $ip]
        }
    }
    return $txt
}


# plot unfolded data and model #
proc az_plot_unfold {args} {
    
    # print help if requested or if no arguments are given #
    if {![llength $args] || [lindex $args 0] == "?" || [lindex $args 0] == "-h"} {
        set help "\n\
    Plot unfolded spectra\n\

    Parameters:\n\
        cmd: either u, eu or eeu, \n\
        out: name of the output file e.g. p for p.plot\n\
        suff: is appended to variable names (_ will be added).\
          default: nothing\n\
        comp: if not 0, plot individual components of the model.\
          default: 0\n\
        resid: if not 0, plot residuals too (resid, del and rat). 
        "
        puts $help
        return
    }

    # check number of arguments #
    if { [llength $args] < 2 || [llength $args] > 5 } {
        puts "\nWrong number of arguments. See -h for help. Quitting ...\n"
        return
    }
    
    chatter 0
    
    # unpack the arguments #
    set cmd  [lindex $args 0]
    set out  [lindex $args 1]
    set suff ""
    if {[llength $args] > 2} {set suff [lindex $args 2]}
    set do_comp 0
    if {[llength $args] > 3} {set do_comp [lindex $args 3]}
    set resid 0
    if {[llength $args] > 4} {set resid 1}
    

    # plot all groups
    set ngrp [tcloutr plotgrp]
    set npar [tcloutr modpar]
    set tmpfile _az_[clock millisec].xcm
    save mod $tmpfile
    
    
    # if doing components; get norms #
    if {$do_comp != 0} {
        set norms {}
        for {set i 1} {$i<=$npar} {incr i} {
            if { [ string equal [lindex [tcloutr pinfo $i] 0] "norm" ] \
              && [ string equal [lindex [tcloutr plink $i] 0] "F" ] \
              && [ expr [lindex [tcloutr par $i] 1] > 0 ] } {
                lappend norms [list $i [lindex [tcloutr param $i] 0]]
            }
        }
    }
    

    set fp [open [format "%s.plot" $out] w ]
    for {set i 1} {$i<=$ngrp} {incr i} {
        puts "@@ plotting data group $i @@"

        set en  [tcloutr plot d x $i]
        set ene [tcloutr plot d xerr $i]
        set d   [tcloutr plot d y $i]
        set de  [tcloutr plot d yerr $i]
        set m   [tcloutr plot d model $i]
        
        
        # residuals if requested #
        if {$resid != 0} {
            set res_1   [tcloutr plot del y $i]
            set res_1_e [tcloutr plot del yerr $i]
            set res_2   [tcloutr plot rat y $i]
            set res_2_e [tcloutr plot rat yerr $i]
            set res_3   [tcloutr plot resid y $i]
            set res_3_e [tcloutr plot resid yerr $i]
        }
        
        
        # get component if requested #
        # first set them = 0; then get one by one #
        set comps {}
        if {$do_comp != 0} {
            foreach n $norms { new [lindex $n 0] 0 .1 -1e20 -1e20 1e20 1e20}
            foreach n $norms {
                new $n
                set mod [tcloutr plot d model $i]
                # this will be 0 if norm doesn't belong to the current ngrp
                if {[::tcl::mathop::+ {*}$mod] != 0 } {
                    lappend comps $mod
                }
                new [lindex $n 0] 0 
            }
        }
        
        
        
        mo po & 0&1
        set u [tcloutr plot $cmd y $i]

        # load model back
        if {[llength $m] > 1} {
            @$tmpfile
        } else {
            mo none
        }

        
        set ss ""
        if {$suff != ""} {set ss [format "%d_%s" $i $suff ]}
        puts $fp "\n\ndescriptor en$ss,+- d$ss,+-"
        for {set ii 0} {$ii<[llength $en]} {incr ii} {
            set fac [expr [lindex $u $ii]/[expr max([lindex $d $ii], 1e-6)] ]
            puts $fp [ format "%3.3e %3.3e %3.3e %3.3e" \
                [lindex $en $ii] [lindex $ene $ii] \
                [lindex $u $ii] [expr [lindex $de $ii]*$fac ] \
            ]
        }
        
        # add model if present #
        if {[llength $m] > 1} {
            puts $fp "\n\ndescriptor m$ss"
            for {set ii 0} {$ii<[llength $en]} {incr ii} {
                set fac [expr [lindex $u $ii]/[expr max([lindex $d $ii], 1e-6)] ]
                puts $fp [ format "%3.3e" [expr [lindex $m $ii]*$fac ]]
            }
        }
        
        # add residuals if present #
        if {$resid != 0} {
            puts $fp "\n\ndescriptor del$ss,+- rat$ss,+- res$ss,+-"
            for {set ii 0} {$ii<[llength $en]} {incr ii} {
                puts $fp [ format "%3.3e 1 %3.3e %3.3e %3.3e %3.3e" [lindex $res_1 $ii] \
                    [lindex $res_2 $ii] [lindex $res_2_e $ii] [lindex $res_3 $ii] \
                    [lindex $res_3_e $ii]]
            }
        }

        
        # Do components if needed #
        if {$do_comp != 0 && [llength $comps] != 0} {
            set desc "\n\ndescriptor "
            for {set ii 1} {$ii<=[llength $comps]} {incr ii} {
                append desc "m${ss}_c$ii "
            }
            puts $fp $desc
            
            for {set ii 0} {$ii<[llength $en]} {incr ii} {
                set fac [expr [lindex $u $ii]/[expr max([lindex $d $ii], 1e-6) ] ]
                set txt ""
                for {set jj 0} {$jj<[llength $comps]} {incr jj} {
                    append txt [ format "%3.3e " \
                        [expr [lindex [lindex $comps $jj] $ii]*$fac ]]
                }
                puts $fp $txt
            }
        }
        
    }
    close $fp
    @$tmpfile
    puts "@@ Done @@"
    chatter 10 10
}


proc q {} {quit}
mdefine log_bpl log(a/e) - log(1 + (e/c)^(-b-1))
mdefine laglogLin a * log(e/eRef) * (f/1e-4)^(-b)
