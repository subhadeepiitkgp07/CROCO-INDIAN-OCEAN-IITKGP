
############ CREATE app.conf for IRENE ###########

if [ ${USE_ATM} -eq 1 ]; then
    echo "${NP_ATM} ./wrfexe" >> app.conf
    if [ ${USE_XIOS_ATM} -eq 1 ]; then
        echo "${NP_XIOS_ATM} ./xios_server.exe" >> app.conf
    fi
fi

if [ ${USE_OCE} -eq 1 ]; then
    if [[ ${MPI_NOLAND} == "TRUE" ]]; then
        echo "${NP_OCE} ./crocox croco.in" >> app.conf
    else
        echo "$(( ${NP_OCEX} * ${NP_OCEY} )) ./crocox croco.in" >> app.conf
    fi
    if [ ${USE_XIOS_OCE} -eq 1 ]; then
        echo "${NP_XIOS_OCE} ./xios_server.exe" >> app.conf
    fi
fi

if [ ${USE_WAV} -eq 1 ]; then
    echo "${NP_WAV} ./wwatch" >> app.conf
fi

if [ ${USE_TOY} -ge 1 ]; then
   for k in `seq 0 $(( ${nbtoy} - 1 ))`; do
       echo "-n ${NP_TOY} ./toy${toytype[$k]}" >> app.conf
   done
fi


