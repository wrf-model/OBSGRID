      program testsndg

c     ... this is a little testing routine that is supposed to generate a 
c         single sounding that the objective analysis program will ingest

c     ... pressure is in Pa, height in m, temperature and dew point are in
c         K, speed is in m/s, and direction is in degrees

c     ... sea level pressure is in Pa, terrain elevation is in m, latitude
c         is in degrees N, longitude is in degrees E

c     ... to put in a surface observation, make only a single level "sounding"
c         and make the value of the height equal the terrain elevation -- PRESTO!

c     ... the first 40 character string may be used for the description of
c         the station (i.e. name city country, etc)

c     ... the second character string we use for our source

c     ... the third string should be left alone, it uses the phrase "FM-35 TEMP"
c         for an upper air station, and should use "FM-12 SYNOP" for surface data

c     ... the fourth string is unused, feel free to experiment with labels!

c     ... bogus data are not subject to quality control

      parameter (kx=10)
      dimension p(kx),z(kx),t(kx),td(kx),spd(kx),dir(kx)
      logical bogus

      data p  /  1000.,   850.,   700.,   500.,   400.,   300.,
     *            250.,   200.,   150.,   100./
      data z  /   100.,  1500.,  3000.,  5500.,  7000.,  9000.,
     *          10500., 12000., 13500., 16000./
      data t  /    14.,     6.,    -4.,   -21.,   -32.,   -45.,
     *            -52.,   -57.,   -57.,   -57./
      data td /    13.,     3.,    -9.,   -28.,   -41.,   -55.,
     *            -62.,   -67.,   -67.,   -67./
      data spd/     1.,     3.,     5.,     7.,     9.,    11.,
     *             13.,    15.,   17.,     19./
      data dir/     0.,    30.,    60.,    90.,   120.,   150.,
     *            180.,   210.,  240.,    270./

      data slp/101325./
      data ter/1./
      data xlat/22./
      data xlon/115./
      data mdate /95073018/
      data bogus /.false./

      do 100 k=1,kx
         p(k)=p(k)*100.
         t(k)=t(k)+273.15
         td(k)=td(k)+273.15
100   continue

      if ( k .eq. 1 ) then
         call write_obs (p,z,t,td,spd,dir, 
     *                    slp, ter, xlat, xlon, mdate, kx, 
     *         '99001  Maybe more site info             ',
     *         'SURFACE DATA FROM ??????????? SOURCE    ',
     *         'FM-12 SYNOP                             ',
     *         '                                        ',
     *         bogus , iseq_num , 2 )
      else
         call write_obs (p,z,t,td,spd,dir, 
     *                    slp, ter, xlat, xlon, mdate, kx, 
     *         '99001  Maybe more site info             ',
     *         'SOUNDINGS FROM ????????? SOURCE         ',
     *         'FM-35 TEMP                              ',
     *         '                                        ',
     *         bogus , iseq_num , 2 )
      endif

      stop 99999
      end

      SUBROUTINE write_obs ( p , z , t , td , spd , dir , 
     *                      slp , ter , xlat , xlon , mdate , kx , 
     * string1 , string2 , string3 , string4 , bogus , iseq_num ,
     * iunit )

      dimension p(kx), z(kx),t(kx),td(kx),spd(kx),dir(kx)

      character *20 date_char
      character *40 string1, string2 , string3 , string4
      CHARACTER *84  rpt_format 
      CHARACTER *22  meas_format 
      CHARACTER *14  end_format
      logical bogus


      rpt_format =  ' ( 2f20.5 , 2a40 , ' 
     *             // ' 2a40 , 1f20.5 , 5i10 , 3L10 , ' 
     *             // ' 2i10 , a20 ,  13( f13.5 , i7 ) ) '
      meas_format =  ' ( 10( f13.5 , i7 ) ) '
      end_format = ' ( 3 ( i7 ) ) ' 

      write (date_char(9:16),fmt='(i8.8)') mdate
      if (mdate/1000000 .GT. 70 ) then
         date_char(7:8)='19'
      else
         date_char(7:8)='20'
      endif
      date_char(17:20)='0000'
      date_char(1:6)='      '

      WRITE ( UNIT = iunit , ERR = 19 , FMT = rpt_format ) 
     *        xlat,xlon, string1 , string2 , 
     *        string3 , string4 , ter, kx*6, 0,0,iseq_num,0, 
     *        .true.,bogus,.false., 
     *         -888888, -888888, date_char , 
     *         slp,0,-888888.,0, -888888.,0, -888888.,0, -888888.,0, 
     *               -888888.,0, 
     *               -888888.,0, -888888.,0, -888888.,0, -888888.,0, 
     *               -888888.,0, 
     *               -888888.,0, -888888.,0
   
      do 100 k = 1 , kx
         WRITE ( UNIT = iunit , ERR = 19 , FMT = meas_format ) 
     *          p(k), 0, z(k),0, t(k),0, td(k),0, 
     *          spd(k),0, dir(k),0, 
     *          -888888.,0, -888888.,0,-888888.,0, -888888.,0
100   continue
      WRITE ( UNIT = iunit , ERR = 19 , FMT = meas_format ) 
     * -777777.,0, -777777.,0,float(kx),0,
     * -888888.,0, -888888.,0, -888888.,0, 
     * -888888.,0, -888888.,0, -888888.,0, 
     * -888888.,0
      WRITE ( UNIT = iunit , ERR = 19 , FMT = end_format )  kx, 0, 0

      return
19    continue
      print *,'troubles writing a sounding'
      stop 19
      END
