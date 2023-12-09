!!!!SJT: edited version of xelbdj2_all.txt for the purposes of compiling routines with other Fortran code
!!!!SJT: removed original driver program
!!!!SJT: removed sample data at end of file
!!!!SJT: created the xelbdj2_all_routines module to facilitate access to the routines (now in the module contains block)
!!!!SJT: explicit declarations for serj and uatan functions were removed because module definition already grants access
!!!!SJT: added use statememt for kind_parameters module (in kind_parameters.f90) in xelbdj2_all_routines module for access to portable kind parameters
!!!!SJT: added implicit none statement in module definition, which extends to all procedures in the contains block
!!!!SJT: modified computations of m from mc to reduce truncation error when mc is close to unity
!!!!SJT: modified upper limits on two do loops in celbdj to avoid an index potentially going out of bounds
!!!!SJT: disabled save and data statements in variable declarations for thread safety
!!!!SJT: added elemental keywords to all procedures (removed write statements that do not occur in practice)
!!!!SJT: variable declarations have been modified to use the kind parameters from the kind_parameters module
!!!!SJT: added intent(in) and intent(out) to argument declarations
!!!!SJT: replaced tab characters with spaces
!!!!SJT: in the code comments, "OG" is short for "Original code" and "SJT" indicates a modification
!!!!SJT: otherwise, the routines are unedited
module xelbdj2_all_routines
 ! module containing subroutines for the calculation of associated Legendre elliptic integrals for standard input ranges
 use kind_parameters
 implicit none
 private
 public :: elbdj2 ! compute the associated Legendre integrals of the first, second, and third kinds for standard input ranges
contains
 !---------------------------------------------------------------------------
 elemental subroutine elbdj2(phi,phic,n,mc_qp,b,d,j) !!SJT
 !subroutine elbdj2(phi,phic,n,mc,b,d,j) !!SJT: OG
 !
 !     Simultaneous computation of associate elliptic integrals
 !     of third kind, B(phi|m), D(phi|m), and J(phi,n|m)
 !
 !     Reference: Fukushima, T, (2012) J. Comp. Appl. Math., 236, 1961-1975
 !       Precise and fast computation of a general incomplete elliptic integral
 !       of third kind by half and double argument transformations
 !
 !     Modified Version by employing "elsbdj" only
 !
 !     Author: Fukushima, T. <Toshio.Fukushima@nao.ac.jp>
 !
 !     Used subprograms: celbd,celbdj,elsbdj,serbd,serj,uatan
 !
 !     Inputs: phi  = argument                0 <= phi  <= PI/2
 !             phic = complementar argument   0 <= phic <= PI/2
 !             n    = characteristic          0 <= n    <= 1
 !             mc   = complementary parameter 0 <= mc   <= 1
 !
 !     Outputs: b, d, j
 !
 !     CAUTION: phi and phic must satisfy condition, phi + phic = PI/2
 !
 real(dp),intent(in) :: phi,phic,n !!SJT
 real(dp),intent(out) :: b,d,j !!SJT
 real(dp) mc,m,nc,h,c,x,d2,z,bc,dc,jc,sz,t,v,t2 !!SJT
 !real(dp) uatan !!SJT: due to module definition, this routine already has access to uatan function 
 real(qp),intent(in) :: mc_qp !!SJT
 !
 if(phi.lt.1.249d0) then   ! Modified: Pascal Leroy's suggestion 2019/08/12
     !call elsbdj(sin(phi),n,mc,b,d,j) !!SJT: OG
     call elsbdj(sin(phi),n,mc_qp,b,d,j) !!SJT
 else
     !m=1.d0-mc !!SJT: original
     m=real(1.0_qp-mc_qp,dp); mc=real(mc_qp,dp) !!SJT
     nc=1.d0-n
     h=n*nc*(n-m)
     c=sin(phic)
     x=c*c
     d2=mc+m*x
     if(x.lt.0.9d0*d2) then
         z=c/sqrt(d2)
         !call elsbdj(z,n,mc,b,d,j) !!SJT: OG
         call elsbdj(z,n,mc_qp,b,d,j) !!SJT
         !call celbdj(nc,mc,bc,dc,jc) !!SJT: original
         call celbdj(nc,mc_qp,bc,dc,jc) !!SJT
                 sz=z*sqrt(1.d0-x)
                 t=sz/nc
                 b=bc-(b-sz)
                 d=dc-(d+sz)
                 j=jc-(j+uatan(t,h))
         else
                 v=mc*(1.d0-x)
                 if(v.lt.x*d2) then
             !call elsbdj(sqrt(1.d0-c*c),n,mc,b,d,j) !!SJT: OG
             call elsbdj(sqrt(1.d0-c*c),n,mc_qp,b,d,j) !!SJT
         else
             t2=(1.d0-x)/d2
             !call elsbdj(sqrt(1.d0-mc*t2),n,mc,b,d,j) !!SJT: OG
             call elsbdj(sqrt(1.d0-mc*t2),n,mc_qp,b,d,j) !!SJT
             !call celbdj(nc,mc,bc,dc,jc) !!SJT: OG
             call celbdj(nc,mc_qp,bc,dc,jc) !!SJT
                     sz=c*sqrt(t2)
                     t=sz/nc
                     b=bc-(b-sz)
                     d=dc-(d+sz)
                     j=jc-(j+uatan(t,h))
                 endif
         endif
 endif
 return
 end
 !---------------------------------------------------------------------------
 elemental subroutine celbd(mc_qp,elb,eld) !!SJT
 !subroutine celbd(mc,elb,eld) !!SJT: OG
 !
 ! Simultaneous computation of associate complete elliptic integrals
 ! of second kind, B(m) and D(m)
 !
 ! Reference: Fukushima, T (2011) Math. Comp., 80, 1725-1743
 !   Precise and fast computation of general complete elliptic integral
 !   of second kind
 !
 ! Author: Fukushima, T. <Toshio.Fukushima@nao.ac.jp>
 !
 real(dp) mc,elk,ele !!SJT
 real(dp),intent(out) :: elb,eld !!SJT
 real(dp) m,mx,kkc,nome,eec,kec
 real(qp),intent(in) :: mc_qp !!SJT
 !
 real(dp) PIQ,PIHALF,PI,PIINV
 parameter (PIQ=0.78539816339744830961566084581988d0)
 parameter (PIHALF=1.5707963267948966192313216916398d0)
 parameter (PI=3.1415926535897932384626433832795d0)
 parameter (PIINV=0.31830988618379067153776752674503d0)
 real(dp) mcold,elbold,eldold
 !save mcold,elbold,eldold !!SJT: disable for thread safety -- variables are initialized at start of each celbd call
 real(dp) Q1,Q2,Q3,Q4,Q5,Q6,Q7,Q8,Q9,Q10,Q11,Q12,Q13,Q14,Q15,Q16
 parameter (Q1=1.d0/16.d0,Q2=1.d0/32.d0,Q3=21.d0/1024.d0)
 parameter (Q4=31.d0/2048.d0,Q5=6257.d0/524288.d0)
 parameter (Q6=10293.d0/1048576.d0,Q7=279025.d0/33554432.d0)
 parameter (Q8=483127.d0/67108864.d0)
 parameter (Q9=435506703.d0/68719476736.d0)
 parameter (Q10=776957575.d0/137438953472.d0)
 parameter (Q11=22417045555.d0/4398046511104.d0)
 parameter (Q12=40784671953.d0/8796093022208.d0)
 parameter (Q13=9569130097211.d0/2251799813685248.d0)
 parameter (Q14=17652604545791.d0/4503599627370496.d0)
 parameter (Q15=523910972020563.d0/144115188075855872.d0)
 parameter (Q16=976501268709949.d0/288230376151711744.d0)
 real(dp) K1,K2,K3,K4,K5,K6,K7
 parameter (K1=1.d0/4.d0)
 parameter (K2=9.d0/64.d0)
 parameter (K3=25.d0/256.d0)
 parameter (K4=1225.d0/16384.d0)
 parameter (K5=3969.d0/65536.d0)
 parameter (K6=53361.d0/1048576.d0)
 parameter (K7=184041.d0/4194304.d0)
 real(dp) B1,B2,B3,B4,B5,B6,B7,B8
 parameter (B1=1.d0/2.d0)
 parameter (B2=1.d0/16.d0)
 parameter (B3=3.d0/128.d0)
 parameter (B4=25.d0/2048.d0)
 parameter (B5=245.d0/32768.d0)
 parameter (B6=1323.d0/262144.d0)
 parameter (B7=7623.d0/2097152.d0)
 parameter (B8=184041.d0/67108864.d0)
 real(dp) D1,D2,D3,D4,D5,D6,D7,D8
 parameter (D1=1.d0/2.d0)
 parameter (D2=3.d0/16.d0)
 parameter (D3=15.d0/128.d0)
 parameter (D4=175.d0/2048.d0)
 parameter (D5=2205.d0/32768.d0)
 parameter (D6=14553.d0/262144.d0)
 parameter (D7=99099.d0/2097152.d0)
 parameter (D8=2760615.d0/67108864.d0)
 real(dp) logq2,dkkc,dddc,ddc,dele,delb,elk1
 !logical first/.TRUE./    !!!SJT: note that logical variable first is always .true. at the start of each routine call
 logical first !!SJT: initialization within declaration not allowed due to thread safety
 first=.true.  !!SJT: initialize outside of declaration for thread safety
 if(first) then
     first=.FALSE.
         mcold=1.d0
         elbold=PIQ
         eldold=PIQ
 endif
 !m=1.d0-mc !!SJT: original
 m=real(1.0_qp-mc_qp,dp); mc=real(mc_qp,dp) !!SJT
 if(abs(mc-mcold).lt.1.11d-16*mc) then
         elb=elbold
         eld=eldold
 elseif(m.lt.1.11d-16) then
     elb=PIQ
         eld=PIQ
 elseif(mc.lt.1.11d-16) then
     elb=1.d0
         eld=0.3862943611198906188344642429164d0-0.5d0*log(mc)
 elseif(mc.lt.0.1d0) then
     nome=mc*(Q1+mc*(Q2+mc*(Q3+mc*(Q4+mc*(Q5+mc*(Q6 &
         +mc*(Q7+mc*(Q8+mc*(Q9+mc*(Q10+mc*(Q11+mc*(Q12 &
         +mc*(Q13+mc*(Q14+mc*(Q15+mc*Q16))))))))))))))) 
     if(mc.lt.0.01d0) then
         dkkc=mc*(K1+mc*(K2+mc*(K3+mc*(K4+mc*(K5+mc*(K6+mc*K7))))))
         dddc=mc*(D1+mc*(D2+mc*(D3+mc*(D4+mc*(D5+mc*(D6+mc*D7))))))
     else
         mx=mc-0.05d0
 !
 ! (K'-1)/(pi/2)
 !
         dkkc=    0.01286425658832983978282698630501405107893d0 &
             +mx*(0.26483429894479586582278131697637750604652d0 &
             +mx*(0.15647573786069663900214275050014481397750d0 &
             +mx*(0.11426146079748350067910196981167739749361d0 &
             +mx*(0.09202724415743445309239690377424239940545d0 &
             +mx*(0.07843218831801764082998285878311322932444d0 &
             +mx*(0.06935260142642158347117402021639363379689d0 &
             +mx*(0.06293203529021269706312943517695310879457d0 &
             +mx*(0.05821227592779397036582491084172892108196d0 &
             +mx*(0.05464909112091564816652510649708377642504d0 &
             +mx*(0.05191068843704411873477650167894906357568d0 &
             +mx*(0.04978344771840508342564702588639140680363d0 &
             +mx*(0.04812375496807025605361215168677991360500d0 &
             ))))))))))))
 !
 ! (K'-E')/(pi/2)
 !
         dddc=    0.02548395442966088473597712420249483947953d0 &
             +mx*(0.51967384324140471318255255900132590084179d0 &
             +mx*(0.20644951110163173131719312525729037023377d0 &
             +mx*(0.13610952125712137420240739057403788152260d0 &
             +mx*(0.10458014040566978574883406877392984277718d0 &
             +mx*(0.08674612915759188982465635633597382093113d0 &
             +mx*(0.07536380269617058326770965489534014190391d0 &
             +mx*(0.06754544594618781950496091910264174396541d0 &
             +mx*(0.06190939688096410201497509102047998554900d0 &
             +mx*(0.05771071515451786553160533778648705873199d0 &
             +mx*(0.05451217098672207169493767625617704078257d0 &
             +mx*(0.05204028407582600387265992107877094920787d0 &
             +mx*(0.05011532514520838441892567405879742720039d0 &
             ))))))))))))
     endif
         kkc=1.d0+dkkc
         logq2=-0.5d0*log(nome)
         elk=kkc*logq2
         dele=-dkkc/kkc+logq2*dddc
         elk1=elk-1.d0
         delb=(dele-mc*elk1)/m
         elb=1.d0+delb
         eld=elk1-delb
 elseif(m.le.0.01d0) then
         elb=PIHALF*(B1+m*(B2+m*(B3+m*(B4+m*(B5+m*(B6+m*(B7+m*B8)))))))
         eld=PIHALF*(D1+m*(D2+m*(D3+m*(D4+m*(D5+m*(D6+m*(D7+m*D8)))))))
 elseif(m.le.0.1d0) then
         mx=0.95d0-mc
         elb=     0.790401413584395132310045630540381158921005d0 &
             +mx*(0.102006266220019154892513446364386528537788d0 &
             +mx*(0.039878395558551460860377468871167215878458d0 &
             +mx*(0.021737136375982167333478696987134316809322d0 &
             +mx*(0.013960979767622057852185340153691548520857d0 &
             +mx*(0.009892518822669142478846083436285145400444d0 &
             +mx*(0.007484612400663335676130416571517444936951d0 &
             +mx*(0.005934625664295473695080715589652011420808d0 &
             +mx*(0.004874249053581664096949448689997843978535d0 &
             +mx*(0.004114606930310886136960940893002069423559d0 &
             +mx*(0.003550452989196176932747744728766021440856d0 &
             +mx*(0.003119229959988474753291950759202798352266d0 &
             )))))))))))
         eld=     0.800602040206397047799296975176819811774784d0 &
             +mx*(0.313994477771767756849615832867393028789057d0 &
             +mx*(0.205913118705551954501930953451976374435088d0 &
             +mx*(0.157744346538923994475225014971416837073598d0 &
             +mx*(0.130595077319933091909091103101366509387938d0 &
             +mx*(0.113308474489758568672985167742047066367053d0 &
             +mx*(0.101454199173630195376251916342483192174927d0 &
             +mx*(0.0929187842072974367037702927967784464949434d0 &
             +mx*(0.0865653801481680871714054745336652101162894d0 &
             +mx*(0.0817279846651030135350056216958053404884715d0 &
             +mx*(0.0779906657291070378163237851392095284454654d0 &
             +mx*(0.075080426851268007156477347905308063808848d0 &
             )))))))))))
 elseif(m.le.0.2d0) then
         mx=0.85d0-mc
         elb=     0.80102406445284489393880821604009991524037d0 &
             +mx*(0.11069534452963401497502459778015097487115d0 &
             +mx*(0.047348746716993717753569559936346358937777d0 &
             +mx*(0.028484367255041422845322166419447281776162d0 &
             +mx*(0.020277811444003597057721308432225505126013d0 &
             +mx*(0.015965005853099119442287313909177068173564d0 &
             +mx*(0.013441320273553634762716845175446390822633d0 &
             +mx*(0.011871565736951439501853534319081030547931d0 &
             +mx*(0.010868363672485520630005005782151743785248d0 &
             +mx*(0.010231587232710564565903812652581252337697d0 &
             +mx*(0.009849585546666211201566987057592610884309d0 &
             +mx*(0.009656606347153765129943681090056980586986d0 &
             )))))))))))
         eld=     0.834232667811735098431315595374145207701720d0 &
             +mx*(0.360495281619098275577215529302260739976126d0 &
             +mx*(0.262379664114505869328637749459234348287432d0 &
             +mx*(0.223723944518094276386520735054801578584350d0 &
             +mx*(0.206447811775681052682922746753795148394463d0 &
             +mx*(0.199809440876486856438050774316751253389944d0 &
             +mx*(0.199667451603795274869211409350873244844882d0 &
             +mx*(0.204157558868236842039815028663379643303565d0 &
             +mx*(0.212387467960572375038025392458549025660994d0 &
             +mx*(0.223948914061499360356873401571821627069173d0 &
             +mx*(0.238708097425597860161720875806632864507536d0 &
             +mx*(0.256707203545463755643710021815937785120030d0 &
             )))))))))))
 elseif(m.le.0.3d0) then
         mx=0.75d0-mc
         elb=     0.81259777291992049322557009456643357559904d0 &
             +mx*(0.12110961794551011284012693733241967660542d0 &
             +mx*(0.057293376831239877456538980381277010644332d0 &
             +mx*(0.038509451602167328057004166642521093142114d0 &
             +mx*(0.030783430301775232744816612250699163538318d0 &
             +mx*(0.027290564934732526869467118496664914274956d0 &
             +mx*(0.025916369289445198731886546557337255438083d0 &
             +mx*(0.025847203343361799141092472018796130324244d0 &
             +mx*(0.026740923539348854616932735567182946385269d0 &
             +mx*(0.028464314554825704963640157657034405579849d0 &
             +mx*(0.030995446237278954096189768338119395563447d0 &
             +mx*(0.034384369179940975864103666824736551261799d0 &
             +mx*(0.038738002072493935952384233588242422046537d0 &
             ))))))))))))
         eld=     0.873152581892675549645633563232643413901757d0 &
             +mx*(0.420622230667770215976919792378536040460605d0 &
             +mx*(0.344231061559450379368201151870166692934830d0 &
             +mx*(0.331133021818721761888662390999045979071436d0 &
             +mx*(0.345277285052808411877017306497954757532251d0 &
             +mx*(0.377945322150393391759797943135325823338761d0 &
             +mx*(0.427378012464553880508348757311170776829930d0 &
             +mx*(0.494671744307822405584118022550673740404732d0 &
             +mx*(0.582685115665646200824237214098764913658889d0 &
             +mx*(0.695799207728083164790111837174250683834359d0 &
             +mx*(0.840018401472533403272555302636558338772258d0 &
             +mx*(1.023268503573606060588689738498395211300483d0 &
             +mx*(1.255859085136282496149035687741403985044122d0 &
             ))))))))))))
 elseif(m.le.0.4d0) then
         mx=0.65d0-mc
         elb=     0.8253235579835158949845697805395190063745d0 &
             +mx*(0.1338621160836877898575391383950840569989d0 &
             +mx*(0.0710112935979886745743770664203746758134d0 &
             +mx*(0.0541784774173873762208472152701393154906d0 &
             +mx*(0.0494517449481029932714386586401273353617d0 &
             +mx*(0.0502221962241074764652127892365024315554d0 &
             +mx*(0.0547429131718303528104722303305931350375d0 &
             +mx*(0.0627462579270016992000788492778894700075d0 &
             +mx*(0.0746698810434768864678760362745179321956d0 &
             +mx*(0.0914808451777334717996463421986810092918d0 &
             +mx*(0.1147050921109978235104185800057554574708d0 &
             +mx*(0.1465711325814398757043492181099197917984d0 &
             +mx*(0.1902571373338462844225085057953823854177d0 &
             ))))))))))))
         eld=     0.9190270392420973478848471774160778462738d0 &
             +mx*(0.5010021592882475139767453081737767171354d0 &
             +mx*(0.4688312705664568629356644841691659415972d0 &
             +mx*(0.5177142277764000147059587510833317474467d0 &
             +mx*(0.6208433913173031070711926900889045286988d0 &
             +mx*(0.7823643937868697229213240489900179142670d0 &
             +mx*(1.0191145350761029126165253557593691585239d0 &
             +mx*(1.3593452027484960522212885423056424704073d0 &
             +mx*(1.8457173023588279422916645725184952058635d0 &
             +mx*(2.5410717031539207287662105618152273788399d0 &
             +mx*(3.5374046552080413366422791595672470037341d0 &
             +mx*(4.9692960029774259303491034652093672488707d0 &
             +mx*(7.0338228700300311264031522795337352226926d0 &
             +mx*(10.020043225034471401553194050933390974016d0 &
             )))))))))))))
 elseif(m.le.0.5d0) then
         mx=0.55d0-mc
         elb=     0.8394795702706129706783934654948360410325d0 &
             +mx*(0.1499164403063963359478614453083470750543d0 &
             +mx*(0.0908319358194288345999005586556105610069d0 &
             +mx*(0.0803470334833417864262134081954987019902d0 &
             +mx*(0.0856384405004704542717663971835424473169d0 &
             +mx*(0.1019547259329903716766105911448528069506d0 &
             +mx*(0.1305748115336160150072309911623351523284d0 &
             +mx*(0.1761050763588499277679704537732929242811d0 &
             +mx*(0.2468351644029554468698889593583314853486d0 &
             +mx*(0.3564244768677188553323196975301769697977d0 &
             +mx*(0.5270025622301027434418321205779314762241d0 &
             +mx*(0.7943896342593047502260866957039427731776d0 &
             +mx*(1.2167625324297180206378753787253096783993d0 &
             ))))))))))))
         eld=     0.9744043665463696730314687662723484085813d0 &
             +mx*(0.6132468053941609101234053415051402349752d0 &
             +mx*(0.6710966695021669963502789954058993004082d0 &
             +mx*(0.8707276201850861403618528872292437242726d0 &
             +mx*(1.2295422312026907609906452348037196571302d0 &
             +mx*(1.8266059675444205694817638548699906990301d0 &
             +mx*(2.8069345309977627400322167438821024032409d0 &
             +mx*(4.4187893290840281339827573139793805587268d0 &
             +mx*(7.0832360574787653249799018590860687062869d0 &
             +mx*(11.515088120557582942290563338274745712174d0 &
             +mx*(18.931511185999274639516732819605594455165d0 &
             +mx*(31.411996938204963878089048091424028309798d0 &
             +mx*(52.520729454575828537934780076768577185134d0 &
             +mx*(88.384854735065298062125622417251073520996d0 &
             +mx*(149.56637449398047835236703116483092644714d0 &
             +mx*(254.31790843104117434615624121937495622372d0 &
             )))))))))))))))
 elseif(m.le.0.6d0) then
         mx=0.45d0-mc
         elb=     0.8554696151564199914087224774321783838373d0 &
             +mx*(0.1708960726897395844132234165994754905373d0 &
             +mx*(0.1213352290269482260207667564010437464156d0 &
             +mx*(0.1282018835749474096272901529341076494573d0 &
             +mx*(0.1646872814515275597348427294090563472179d0 &
             +mx*(0.2374189087493817423375114793658754489958d0 &
             +mx*(0.3692081047164954516884561039890508294508d0 &
             +mx*(0.6056587338479277173311618664015401963868d0 &
             +mx*(1.0337055615578127436826717513776452476106d0 &
             +mx*(1.8189884893632678849599091011718520567105d0 &
             +mx*(3.2793776512738509375806561547016925831128d0 &
             +mx*(6.0298883807175363312261449542978750456611d0 &
             +mx*(11.269796855577941715109155203721740735793d0 &
             +mx*(21.354577850382834496786315532111529462693d0 &
             )))))))))))))
         eld=     1.04345529511513353426326823569160142342838d0 &
             +mx*(0.77962572192850485048535711388072271372632d0 &
             +mx*(1.02974236093206758187389128668777397528702d0 &
             +mx*(1.62203722341135313022433907993860147395972d0 &
             +mx*(2.78798953118534762046989770119382209443756d0 &
             +mx*(5.04838148737206914685643655935236541332892d0 &
             +mx*(9.46327761194348429539987572314952029503864d0 &
             +mx*(18.1814899494276679043749394081463811247757d0 &
             +mx*(35.5809805911791687037085198750213045708148d0 &
             +mx*(70.6339354619144501276254906239838074917358d0 &
             +mx*(141.828580083433059297030133195739832297859d0 &
             +mx*(287.448751250132166257642182637978103762677d0 &
             +mx*(587.115384649923076181773192202238389711345d0 &
             +mx*(1207.06543522548061603657141890778290399603d0 &
             +mx*(2495.58872724866422273012188618178997342537d0 &
             +mx*(5184.69242939480644062471334944523925163600d0 &
             +mx*(10817.2133369041327524988910635205356016939d0 &
             ))))))))))))))))
 elseif(m.le.0.7d0) then
         mx=0.35d0-mc
         elb=     0.8739200618486431359820482173294324246058d0 &
             +mx*(0.1998140574823769459497418213885348159654d0 &
             +mx*(0.1727696158780152128147094051876565603862d0 &
             +mx*(0.2281069132842021671319791750725846795701d0 &
             +mx*(0.3704681411180712197627619157146806221767d0 &
             +mx*(0.6792712528848205545443855883980014994450d0 &
             +mx*(1.3480084966817573020596179874311042267679d0 &
             +mx*(2.8276709768538207038046918250872679902352d0 &
             +mx*(6.1794682501239140840906583219887062092430d0 &
             +mx*(13.935686010342811497608625663457407447757d0 &
             +mx*(32.218929281059722026322932181420383764028d0 &
             +mx*(76.006962959226101026399085304912635262362d0 &
             +mx*(182.32144908775406957609058046006949657416d0 &
             +mx*(443.51507644112648158679360783118806161062d0 &
             +mx*(1091.8547229028388292980623647414961662223d0 &
             +mx*(2715.7658664038195881056269799613407111521d0 &
             )))))))))))))))
         eld=     1.13367833657573316571671258513452768536080d0 &
             +mx*(1.04864317372997039116746991765351150490010d0 &
             +mx*(1.75346504119846451588826580872136305225406d0 &
             +mx*(3.52318272680338551269021618722443199230946d0 &
             +mx*(7.74947641381397458240336356601913534598302d0 &
             +mx*(17.9864500558507330560532617743406294626849d0 &
             +mx*(43.2559163462326133313977294448984936591235d0 &
             +mx*(106.681534454096017031613223924991564429656d0 &
             +mx*(268.098486573117433951562111736259672695883d0 &
             +mx*(683.624114850289804796762005964155730439745d0 &
             +mx*(1763.49708521918740723028849567007874329637d0 &
             +mx*(4592.37475383116380899419201719007475759114d0 &
             +mx*(12053.4410190488892782190764838488156555734d0 &
             +mx*(31846.6630207420816960681624497373078887317d0 &
             +mx*(84621.2213590568080177035346867495326879117d0 &
             +mx*(225956.423182907889987641304430180593010940d0 &
             +mx*(605941.517281758859958050194535269219533685d0 &
             +mx*(1.63108259953926832083633544697688841456604d6 &
             )))))))))))))))))
 elseif(m.le.0.8d0) then
         mx=0.25d0-mc
         elb=     0.895902820924731621258525533131864225704d0 &
             +mx*(0.243140003766786661947749288357729051637d0 &
             +mx*(0.273081875594105531575351304277604081620d0 &
             +mx*(0.486280007533573323895498576715458103274d0 &
             +mx*(1.082747437228230914750752674136983406683d0 &
             +mx*(2.743445290986452500459431536349945437824d0 &
             +mx*(7.555817828670234627048618342026400847824d0 &
             +mx*(22.05194082493752427472777448620986154515d0 &
             +mx*(67.15640644740229407624192175802742979626d0 &
             +mx*(211.2722537881770961691291434845898538537d0 &
             +mx*(681.9037843053270682273212958093073895805d0 &
             +mx*(2246.956231592536516768812462150619631201d0 &
             +mx*(7531.483865999711792004783423815426725079d0 &
             +mx*(25608.51260130241579018675054866136922157d0 &
             +mx*(88140.74740089604971425934283371277143256d0 &
             +mx*(306564.4242098446591430938434419151070722d0 &
             +mx*(1.076036077811072193752770590363885180738d6 &
             +mx*(3.807218502573632648224286313875985190526d6 &
             +mx*(1.356638224422139551020110323739879481042d7 &
             ))))))))))))))))))
         eld=     1.26061282657491161418014946566845780315983d0 &
             +mx*(1.54866563808267658056930177790599939977154d0 &
             +mx*(3.55366941187160761540650011660758187283401d0 &
             +mx*(9.90044467610439875577300608183010716301714d0 &
             +mx*(30.3205666174524719862025105898574414438275d0 &
             +mx*(98.1802586588830891484913119780870074464833d0 &
             +mx*(329.771010434557055036273670551546757245808d0 &
             +mx*(1136.65598974289039303581967838947708073239d0 &
             +mx*(3993.83433574622979757935610692842933356144d0 &
             +mx*(14242.7295865552708506232731633468180669284d0 &
             +mx*(51394.7572916887209594591528374806790960057d0 &
             +mx*(187246.702914623152141768788258141932569037d0 &
             +mx*(687653.092375389902708761221294282367947659d0 &
             +mx*(2.54238553565398227033448846432182516906624d6 &
             +mx*(9.45378121934749027243313241962076028066811d6 &
             +mx*(3.53283630179709170835024033154326126569613d7 &
             +mx*(1.32593262383393014923560730485845833322771d8 &
             +mx*(4.99544968184054821463279808395426941549833d8 &
             +mx*(1.88840934729443872364972817525484292678543d9 &
             +mx*(7.16026753447893719179055010636502508063102d9 &
             +mx*(2.72233079469633962247554894093591262281929d10 &
         ))))))))))))))))))))
 elseif(m.le.0.85d0) then
         mx=0.175d0-mc
         elb=     0.915922052601931494319853880201442948834592d0 &
             +mx*(0.294714252429483394379515488141632749820347d0 &
             +mx*(0.435776709264636140422971598963772380161131d0 &
             +mx*(1.067328246493644238508159085364429570207744d0 &
             +mx*(3.327844118563268085074646976514979307993733d0 &
             +mx*(11.90406004445092906188837729711173326621810d0 &
             +mx*(46.47838820224626393512400481776284680677096d0 &
             +mx*(192.7556002578809476962739389101964074608802d0 &
             +mx*(835.3356299261900063712302517586717381557137d0 &
             +mx*(3743.124548343029102644419963712353854902019d0 &
             +mx*(17219.07731004063094108708549153310467326395d0 &
             +mx*(80904.60401669850158353080543152212152282878d0 &
             +mx*(386808.3292751742460123683674607895217760313d0 &
             +mx*(1.876487670110449342170327796786290400635732d6 &
             +mx*(9.216559908641567755240142886998737950775910d6 &
             ))))))))))))))
         eld=     1.402200569110579095046054435635136986038164d0 &
             +mx*(2.322205897861749446534352741005347103992773d0 &
             +mx*(7.462158366466719682730245467372788273333992d0 &
             +mx*(29.43506890797307903104978364254987042421285d0 &
             +mx*(128.1590924337895775262509354898066132182429d0 &
             +mx*(591.0807036911982326384997979640812493154316d0 &
             +mx*(2830.546229607726377048576057043685514661188d0 &
             +mx*(13917.76431889392229954434840686741305556862d0 &
             +mx*(69786.10525163921228258055074102587429394212d0 &
             +mx*(355234.1420341879634781808899208309503519936d0 &
             +mx*(1.830019186413931053503912913904321703777885d6 &
             +mx*(9.519610812032515607466102200648641326190483d6 &
             +mx*(4.992086875574849453986274042758566713803723d7 &
             +mx*(2.635677009826023473846461512029006874800883d8 &
             +mx*(1.399645765120061118824228996253541612110338d9 &
             +mx*(7.469935792837635004663183580452618726280406d9 &
             +mx*(4.004155595835610574316003488168804738481448d10 &
             +mx*(2.154630668144966654449602981243932210695662d11 &
             )))))))))))))))))
 else
     mx=0.125d0-mc
         elb=     0.931906061029524827613331428871579482766771d0 &
             +mx*(0.348448029538453860999386797137074571589376d0 &
             +mx*(0.666809178846938247558793864839434184202736d0 &
             +mx*(2.210769135708128662563678717558631573758222d0 &
             +mx*(9.491765048913406881414290930355300611703187d0 &
             +mx*(47.09304791027740853381457907791343619298913d0 &
             +mx*(255.9200460211233087050940506395442544885608d0 &
             +mx*(1480.029532675805407554800779436693505109703d0 &
             +mx*(8954.040904734313578374783155553041065984547d0 &
             +mx*(56052.48220982686949967604699243627330816542d0 &
             +mx*(360395.7241626000916973524840479780937869149d0 &
             +mx*(2.367539415273216077520928806581689330885103d6 &
             +mx*(1.582994957277684102454906900025484391190264d7 &
             +mx*(1.074158093278511100137056972128875270067228d8 &
             +mx*(7.380585460239595691878086073095523043390649d8 &
             +mx*(5.126022002555101496684687154904781856830296d9 &
             +mx*(3.593534065502416588712409180013118409428367d10 &
             +mx*(2.539881257612812212004146637239987308133582d11 &
             +mx*(1.808180007145359569674767150594344316702507d12 &
         ))))))))))))))))))
         eld=     1.541690112721819084362258323861459983048179d0 &
             +mx*(3.379176214579645449453938918349243359477706d0 &
             +mx*(14.94058385670236671625328259137998668324435d0 &
             +mx*(81.91773929235074880784578753539752529822986d0 &
             +mx*(497.4900546551479866036061853049402721939835d0 &
             +mx*(3205.184010234846235275447901572262470252768d0 &
             +mx*(21457.32237355321925571253220641357074594515d0 &
             +mx*(147557.0156564174712105689758692929775004292d0 &
             +mx*(1.035045290185256525452269053775273002725343d6 &
             +mx*(7.371922334832212125197513363695905834126154d6 &
             +mx*(5.314344395142401141792228169170505958906345d7 &
             +mx*(3.868823475795976312985118115567305767603128d8 &
             +mx*(2.839458401528033778425531336599562337200510d9 &
             +mx*(2.098266122943898941547136470383199468548861d10 &
             +mx*(1.559617754017662417944194874282275405676282d11 &
             +mx*(1.165096220419884791236699872205721392201682d12 &
             +mx*(8.742012983013913804987431275193291316808766d12 &
             +mx*(6.584725462672366918676967847406180155459650d13 &
             +mx*(4.976798737062434393396993620379481464465749d14 &
             +mx*(3.773018634056605404718444239040628892506293d15 &
             +mx*(2.868263194837819660109735981973458220407767d16 &
         ))))))))))))))))))))
 endif
 mcold=mc
 elbold=elb
 eldold=eld
 return
 end
 !---------------------------------------------------------------------------
 elemental subroutine celbdj(nc0,mc0_qp,celb,celd,celj)
 !subroutine celbdj(nc0,mc0,celb,celd,celj) !!SJT: OG
 !
 ! Simultaneous computation of associate complete elliptic integrals
 ! of third kind, B(m), D(m), and J(n|m)
 !
 ! Reference: Fukushima, T, (2013) J. Comp. Appl. Math., 253, 142-157
 !     Fast computation of a general complete elliptic integral
 !     of third kind by half and double argument transformations
 !
 ! Restricted Version with Assumptions: 0 < mc, 0 < nc
 !
 ! Author: Fukushima, T. <Toshio.Fukushima@nao.ac.jp>
 !
 real(dp),intent(in) :: nc0
 real(dp),intent(out) :: celb,celd,celj !!SJT
 logical flag
 integer(isp) IMAX,i,is,ie,err
 parameter (IMAX=40)
 real(dp) y(0:IMAX),x(0:IMAX),c(0:IMAX),d(0:IMAX),a(0:IMAX)
 real(dp) mc,mc0,nc,m,n,yj,celk,yi,ye,dj,m1,kc0,temp !!SJT
 real(qp),intent(in) :: mc0_qp !!SJT
 real(qp) :: mc_qp !!SJT
 real(dp) PIHALF,EPS,THIRD
 parameter (PIHALF=1.5707963267948966d0)
 parameter (EPS=1.11d-16)
 parameter (THIRD=1.d0/3.d0)
 real(dp) B(IMAX)
 !logical first /.TRUE./  !!!SJT: note that logical variable first is .true. at the start of each routine call
 !save B,first !!SJT: disable for thread safety
 logical first !!SJT: initialization cannot be in declaration for thread-safe calculations
 first=.true. !!SJT: initialization outside of declaration
 mc0=real(mc0_qp,dp) !!SJT
 !if(mc0.le.0.d0) then !!SJT: does not occur in practice
 !    write(*,*) "(celbdj) Out of domain: mc <= 0"
 !    return
 !endif
 !if(nc0.le.0.d0) then
 !    write(*,*) "(celbdj) Out of domain: nc <= 0"
 !    return
 !endif
 if(mc0.lt.1.d0) then
     !mc=mc0;nc=nc0
     mc_qp=mc0_qp; nc=nc0 !!SJT
 elseif(mc0.gt.1.d0) then
     !mc=1.d0/mc0;nc=nc0*mc !!SJT: OG
     mc_qp=1.0_qp/mc0_qp; nc=nc0*real(mc_qp,dp) !!SJT
 else
     celb=PIHALF*0.5d0
     celd=PIHALF*0.5d0
     celj=PIHALF/(nc0+sqrt(nc0))
     return
 endif
 if(first) then
     first=.FALSE.
     do i=1,IMAX
         B(i)=1.d0/dble(2*i+1)
     enddo
 endif
 mc=real(mc_qp,dp) !!SJT
 !call celbd(mc,celb,celd) !!SJT: OG
 call celbd(mc_qp,celb,celd) !!SJT
 if(nc.eq.mc) then
     celj=celb/mc
     goto 4  
 endif
 !m=1.d0-mc; n=1.d0-nc !!SJT: original
 m=real(1.0_qp-mc_qp,dp); n=1.d0-nc !!SJT
 flag=nc.lt.mc.or.(n*nc).gt.(nc-mc)
 if(flag) then
     y(0)=(nc-mc)/(nc*m)
 else
     y(0)=n/m
 endif
 is=0
 if(y(0).gt.0.5d0) then
     x(0)=1.d0-y(0)
     do i=0,IMAX-1 !!SJT
     !do i=0,IMAX !!SJT: OG
         c(i)=sqrt(x(i))
         d(i)=sqrt(mc+m*x(i))
         x(i+1)=(c(i)+d(i))/(1.d0+d(i))
         if(x(i+1).gt.0.5d0) then
             y(i+1)=y(i)/((1.d0+c(i))*(1.d0+d(i)))
             is=i+1
             goto 1
         endif
         y(i+1)=1.d0-x(i+1)
     enddo
     !write(*,"(a30,i5)") "(celbdj) No Conv. x-Transf. i=",i !!SJT: does not occur in practice
     return
 endif
 1 continue
 do i=is,IMAX-1 !!SJT
 !do i=is,IMAX !!SJT: OG
     c(i)=sqrt(1.d0-y(i))
     d(i)=sqrt(1.d0-m*y(i))
     y(i+1)=y(i)/((1.d0+c(i))*(1.d0+d(i)))
     if(abs(y(i+1)).lt.0.325d0) goto 2
 enddo
 !write(*,"(a30,i5)") "(celbdj) No Conv. y-Transf. i=",i !!SJT: does not occur in practice
 return
 2 continue
 ie=i+1
 ye=y(ie)
 celk=celb+celd
 a(0)=celd
 celj=a(0)
 yi=ye
 a(1)=((1.d0+2.d0*m)*celd-celb)*THIRD
 dj=a(1)*yi
 i=1
 if(abs(dj).lt.EPS*abs(celj)) goto 3
 celj=celj+dj
 m1=1.d0+m
 do i=2,IMAX
     yi=yi*ye
     a(i)=(1.d0-B(i))*m1*a(i-1)-(1.d0-2.d0*B(i))*m*a(i-2)
     dj=a(i)*yi
     if(abs(dj).lt.EPS*abs(celj)) goto 3
     celj=celj+dj
 enddo
 !write(*,"(a30,i5)") "(celbdj) No Conv. Series. i=",i !!SJT: does not occur in practice
 return
 3 continue
 do i=ie-1,0,-1
     celj=(2.d0*(c(i)+d(i))*celj-y(i)*celk)/(c(i)*d(i)*(1.d0+c(i))*(1.d0+d(i)))
 enddo
 if(flag) then
     celj=(nc*celk-mc*celj)/(nc*nc)
 endif
 4 continue
 if(mc0.gt.1.d0) then
     kc0=sqrt(mc0)
     temp=celb
     celb=celd/kc0
     celd=temp/kc0
     celj=celj/(mc0*kc0)
 endif
 !write(*,"(a30,1p3e15.7)") "(celbdj) nc0,mc0,celj=",nc0,mc0,celj
 return
 end
 !---------------------------------------------------------------------------
 elemental subroutine elsbdj(s0,n,mc_qp,b,d,j) !!SJT
 !subroutine elsbdj(s0,n,mc,b,d,j) !!SJT: OG
 !
 ! Simultaneous computation of associate elliptic integrals
 ! of third kind, B(phi|m), D(phi|m), and J(phi,n|m)
 ! by using the half/double argument transformation of sn functions
 !
 ! Reference: Fukushima, T, (2012) J. Comp. Appl. Math., 236, 1961-1975
 !   Precise and fast computation of a general incomplete elliptic integral
 !   of third kind by half and double argument transformations
 !
 real(dp),intent(in) :: s0,n !!SJT
 real(dp),intent(out) :: b,d,j !!SJT
 real(dp) m,mc,h,del,s,y,c,sy,t !!SJT
 real(qp),intent(in) :: mc_qp !!SJT
 real(dp) yy(11),ss(11),cd(11)
 !real(dp) serj,uatan !!SJT: due to module definition, this routine already has access to serj and uatan functions
 integer(isp) i,k
 
 ! write(*,*) "(elsj) s0,n,mc=",s0,n,mc
 
 !m=1.d0-mc !!SJT: original
 m=real(1.0_qp-mc_qp,dp); mc=real(mc_qp,dp) !!SJT
 !kc=mc**0.5d0; m=(1.d0+kc)*(1.d0-kc) !!SJT: to reduce truncation error when mc is close to unity
 h=n*(1.d0-n)*(n-m)
 !del=0.04094d0-0.00652d0*m
 !del=0.03429d0
 !del=0.02448d0  ! JA
 del=0.01622d0   ! J9
 !del=0.00969d0  ! J8
 !del=0.00500d0  ! J7
 !del=2.073d-3    ! J6
 !del=6.037d-4    ! J5
 
 s=s0
 y=s*s
 if(y.lt.del) then
         call serbd(y,m,b,d)
         b=s*b
         d=s*y*d
         j=s*serj(y,n,m)
 ! write(*,"(a20,1p3e10.2)") "(elsbdj) b,d,j=",b,d,j
         return
 endif
 yy(1)=y
 ss(1)=s
 ! write(*,"(a20,i10,1pe10.2)") "(elsbdj) i,y=",1,y
 do i=1,10
     c=sqrt(1.d0-y)
     d=sqrt(1.d0-m*y)
     y=y/((1.d0+c)*(1.d0+d))
         yy(i+1)=y
         ss(i+1)=sqrt(y)
         cd(i)=c*d
 ! write(*,"(a20,i10,1pe10.2)") "(elsbdj) i,y=",i+1,y
         if(y.lt.del) then
                 goto 1
         endif
 enddo
 !write(*,*) "(elsbdj) too many iterations: s0,m=",s0,m SJT: does not occur in practice
 1 continue
 ! write(*,"(a20,i10,1pe10.2)") "(elsbdj) i,y=",i+1,y
 call serbd(y,m,b,d)
 b=ss(i+1)*b
 d=ss(i+1)*y*d
 j=ss(i+1)*serj(y,n,m)
 do k=i,1,-1
         sy=ss(k)*yy(k+1)
         t=sy/(1.d0-n*(yy(k)-yy(k+1)*cd(k)))
         b=2.d0*b-sy
         d=d+(d+sy)
         j=j+(j+uatan(t,h))
 enddo
 ! write(*,"(a20,1p3e10.2)") "(elsbdj) b,d,j=",b,d,j
 return
 end
 !---------------------------------------------------------------------------
 elemental subroutine serbd(y,m,b,d)
 !
 ! Simultaneous computation of associate elliptic integrals,
 ! B(phi|m) and D(phi|m), for small arguments by the series expansion
 !
 ! Reference: Fukushima, T (2012) J. Comp. Appl. Math., 235, 4140-4148
 !   Precise and fast computation of general incomplete elliptic integral
 !   of second kind by half and double argument transformations
 !
 real(dp),intent(in) :: y,m !!SJT
 real(dp),intent(out) :: b,d !!SJT
 !
 real(dp) F1,F2,F3,F4
 real(dp) F10,F20,F21,F30,F31,F40,F41,F42
 real(dp) F5,F50,F51,F52,F6,F60,F61,F62,F63
 real(dp) F7,F70,F71,F72,F73,F8,F80,F81,F82,F83,F84
 real(dp) F9,F90,F91,F92,F93,F94
 real(dp) FA,FA0,FA1,FA2,FA3,FA4,FA5
 real(dp) FB,FB0,FB1,FB2,FB3,FB4,FB5
 parameter (F10=1.d0/6.d0)
 parameter (F20=3.d0/40.d0)
 parameter (F21=2.d0/40.d0)
 parameter (F30=5.d0/112.d0)
 parameter (F31=3.d0/112.d0)
 parameter (F40=35.d0/1152.d0)
 parameter (F41=20.d0/1152.d0)
 parameter (F42=18.d0/1152.d0)
 parameter (F50=63.d0/2816.d0)
 parameter (F51=35.d0/2816.d0)
 parameter (F52=30.d0/2816.d0)
 parameter (F60=231.d0/13312.d0)
 parameter (F61=126.d0/13312.d0)
 parameter (F62=105.d0/13312.d0)
 parameter (F63=100.d0/13312.d0)
 parameter (F70=429.d0/30720.d0)
 parameter (F71=231.d0/30720.d0)
 parameter (F72=189.d0/30720.d0)
 parameter (F73=175.d0/30720.d0)
 parameter (F80=6435.d0/557056.d0)
 parameter (F81=3432.d0/557056.d0)
 parameter (F82=2722.d0/557056.d0)
 parameter (F83=2520.d0/557056.d0)
 parameter (F84=2450.d0/557056.d0)
 parameter (F90=12155.d0/1245184.d0)
 parameter (F91=6435.d0/1245184.d0)
 parameter (F92=5148.d0/1245184.d0)
 parameter (F93=4620.d0/1245184.d0)
 parameter (F94=4410.d0/1245184.d0)
 parameter (FA0=46189.d0/5505024.d0)
 parameter (FA1=24310.d0/5505024.d0)
 parameter (FA2=19305.d0/5505024.d0)
 parameter (FA3=17160.d0/5505024.d0)
 parameter (FA4=16170.d0/5505024.d0)
 parameter (FA5=15876.d0/5505024.d0)
 parameter (FB0=88179.d0/12058624.d0)
 parameter (FB1=46189.d0/12058624.d0)
 parameter (FB2=36465.d0/12058624.d0)
 parameter (FB3=32175.d0/12058624.d0)
 parameter (FB4=30030.d0/12058624.d0)
 parameter (FB5=29106.d0/12058624.d0)
 real(dp) A1,A2,A3,A4,A5,A6,A7,A8,A9,AA,AB
 parameter (A1=3.d0/5.d0)
 parameter (A2=5.d0/7.d0)
 parameter (A3=7.d0/9.d0)
 parameter (A4=9.d0/11.d0)
 parameter (A5=11.d0/13.d0)
 parameter (A6=13.d0/15.d0)
 parameter (A7=15.d0/17.d0)
 parameter (A8=17.d0/19.d0)
 parameter (A9=19.d0/21.d0)
 parameter (AA=21.d0/23.d0)
 parameter (AB=23.d0/25.d0)
 real(dp) B1,B2,B3,B4,B5,B6,B7,B8,B9,BA,BB
 real(dp) D0,D1,D2,D3,D4,D5,D6,D7,D8,D9,DA,DB
 parameter (D0=1.d0/3.d0)
 !
 !        write(*,*) "(serbd) y,m=",y,m
 F1=F10+m*F10
 F2=F20+m*(F21+m*F20)
 F3=F30+m*(F31+m*(F31+m*F30))
 F4=F40+m*(F41+m*(F42+m*(F41+m*F40)))
 F5=F50+m*(F51+m*(F52+m*(F52+m*(F51+m*F50))))
 F6=F60+m*(F61+m*(F62+m*(F63+m*(F62+m*(F61+m*F60)))))
 F7=F70+m*(F71+m*(F72+m*(F73+m*(F73+m*(F72+m*(F71+m*F70))))))
 F8=F80+m*(F81+m*(F82+m*(F83+m*(F84+m*(F83+m*(F82+m*(F81+m*F80)))))))
 F9=F90+m*(F91+m*(F92+m*(F93+m*(F94+m*(F94+m*(F93 &
     +m*(F92+m*(F91+m*F90))))))))
 FA=FA0+m*(FA1+m*(FA2+m*(FA3+m*(FA4+m*(FA5+m*(FA4 &
     +m*(FA3+m*(FA2+m*(FA1+m*FA0)))))))))
 FB=FB0+m*(FB1+m*(FB2+m*(FB3+m*(FB4+m*(FB5+m*(FB5+ &
     m*(FB4+m*(FB3+m*(FB2+m*(FB1+m*FB0))))))))))
 !
 D1=F1*A1
 D2=F2*A2
 D3=F3*A3
 D4=F4*A4
 D5=F5*A5
 D6=F6*A6
 D7=F7*A7
 D8=F8*A8
 D9=F9*A9
 DA=FA*AA
 DB=FB*AB
 d=D0+y*(D1+y*(D2+y*(D3+y*(D4+y*(D5+y*(D6+y*(D7+y*(D8 &
     +y*(D9+y*(DA+y*DB))))))))))
 B1=F1-D0
 B2=F2-D1
 B3=F3-D2
 B4=F4-D3
 B5=F5-D4
 B6=F6-D5
 B7=F7-D6
 B8=F8-D7
 B9=F9-D8
 BA=FA-D9
 BB=FB-DA
 b=1.d0+y*(B1+y*(B2+y*(B3+y*(B4+y*(B5+y*(B6+y*(B7+y*(B8 &
     +y*(B9+y*(BA+y*BB))))))))))
 return
 end
 !---------------------------------------------------------------------------
 elemental real(dp) function serj(y,n,m)
 !
 ! Computation of associate elliptic integral J(phi,n|m)
 ! for small arguments by the series expansion
 !
 ! Reference: Fukushima, T, (2012) J. Comp. Appl. Math., 236, 1961-1975
 !   Precise and fast computation of a general incomplete elliptic integral
 !   of third kind by half and double argument transformations
 !
 real(dp),intent(in) :: y,n,m !!SJT
 
 real(dp) J1,J2,J3,J4,J5,J6,J7,J8,J9,JA
 
 real(dp) J100,J200,J201,J210,J300,J301,J302,J310,J311,J320
 real(dp) J400,J401,J402,J403,J410,J411,J412,J420,J421,J430
 real(dp) J500,J501,J502,J503,J504,J510,J511,J512,J513,J520
 real(dp) J521,J522,J530,J531,J540
 real(dp) J600,J601,J602,J603,J604,J605,J610,J611,J612,J613,J614
 real(dp) J620,J621,J622,J623,J630,J631,J632,J640,J641,J650
 real(dp) J700,J701,J702,J703,J704,J705,J706
 real(dp) J710,J711,J712,J713,J714,J715,J720,J721,J722,J723,J724
 real(dp) J730,J731,J732,J733,J740,J741,J742,J750,J751,J760
 real(dp) J800,J801,J802,J803,J804,J805,J806,J807
 real(dp) J810,J811,J812,J813,J814,J815,J816
 real(dp) J820,J821,J822,J823,J824,J825,J830,J831,J832,J833,J834
 real(dp) J840,J841,J842,J843,J850,J851,J852,J860,J861,J870
 real(dp) J900,J901,J902,J903,J904,J905,J906,J907,J908
 real(dp) J910,J911,J912,J913,J914,J915,J916,J917
 real(dp) J920,J921,J922,J923,J924,J925,J926
 real(dp) J930,J931,J932,J933,J934,J935,J940,J941,J942,J943,J944
 real(dp) J950,J951,J952,J953,J960,J961,J962,J970,J971,J980
 real(dp) JA00,JA01,JA02,JA03,JA04,JA05,JA06,JA07,JA08,JA09
 real(dp) JA10,JA11,JA12,JA13,JA14,JA15,JA16,JA17,JA18
 real(dp) JA20,JA21,JA22,JA23,JA24,JA25,JA26,JA27
 real(dp) JA30,JA31,JA32,JA33,JA34,JA35,JA36
 real(dp) JA40,JA41,JA42,JA43,JA44,JA45,JA50,JA51,JA52,JA53,JA54
 real(dp) JA60,JA61,JA62,JA63,JA70,JA71,JA72,JA80,JA81,JA90
 
 parameter (J100=1.d0/3.d0)
 
 parameter (J200=1.d0/10.d0)
 parameter (J201=2.d0/10.d0)
 parameter (J210=1.d0/10.d0)
 
 parameter (J300=3.d0/56.d0)
 parameter (J301=4.d0/56.d0)
 parameter (J302=8.d0/56.d0)
 parameter (J310=2.d0/56.d0)
 parameter (J311=4.d0/56.d0)
 parameter (J320=3.d0/56.d0)
 
 parameter (J400=5.d0/144.d0)
 parameter (J401=6.d0/144.d0)
 parameter (J402=8.d0/144.d0)
 parameter (J403=16.d0/144.d0)
 parameter (J410=3.d0/144.d0)
 parameter (J411=4.d0/144.d0)
 parameter (J412=8.d0/144.d0)
 parameter (J420=3.d0/144.d0)
 parameter (J421=6.d0/144.d0)
 parameter (J430=5.d0/144.d0)
 
 parameter (J500=35.d0/1408.d0)
 parameter (J501=40.d0/1408.d0)
 parameter (J502=48.d0/1408.d0)
 parameter (J503=64.d0/1408.d0)
 parameter (J504=128.d0/1408.d0)
 parameter (J510=20.d0/1408.d0)
 parameter (J511=24.d0/1408.d0)
 parameter (J512=32.d0/1408.d0)
 parameter (J513=64.d0/1408.d0)
 parameter (J520=18.d0/1408.d0)
 parameter (J521=24.d0/1408.d0)
 parameter (J522=48.d0/1408.d0)
 parameter (J530=20.d0/1408.d0)
 parameter (J531=40.d0/1408.d0)
 parameter (J540=35.d0/1408.d0)
 
 parameter (J600=63.d0/3328.d0)
 parameter (J601=70.d0/3328.d0)
 parameter (J602=80.d0/3328.d0)
 parameter (J603=96.d0/3328.d0)
 parameter (J604=128.d0/3328.d0)
 parameter (J605=256.d0/3328.d0)
 parameter (J610=35.d0/3328.d0)
 parameter (J611=40.d0/3328.d0)
 parameter (J612=48.d0/3328.d0)
 parameter (J613=64.d0/3328.d0)
 parameter (J614=128.d0/3328.d0)
 parameter (J620=30.d0/3328.d0)
 parameter (J621=36.d0/3328.d0)
 parameter (J622=48.d0/3328.d0)
 parameter (J623=96.d0/3328.d0)
 parameter (J630=30.d0/3328.d0)
 parameter (J631=40.d0/3328.d0)
 parameter (J632=80.d0/3328.d0)
 parameter (J640=35.d0/3328.d0)
 parameter (J641=70.d0/3328.d0)
 parameter (J650=63.d0/3328.d0)
 
 parameter (J700=231.d0/15360.d0)
 parameter (J701=252.d0/15360.d0)
 parameter (J702=280.d0/15360.d0)
 parameter (J703=320.d0/15360.d0)
 parameter (J704=384.d0/15360.d0)
 parameter (J705=512.d0/15360.d0)
 parameter (J706=1024.d0/15360.d0)
 parameter (J710=126.d0/15360.d0)
 parameter (J711=140.d0/15360.d0)
 parameter (J712=160.d0/15360.d0)
 parameter (J713=192.d0/15360.d0)
 parameter (J714=256.d0/15360.d0)
 parameter (J715=512.d0/15360.d0)
 parameter (J720=105.d0/15360.d0)
 parameter (J721=120.d0/15360.d0)
 parameter (J722=144.d0/15360.d0)
 parameter (J723=192.d0/15360.d0)
 parameter (J724=384.d0/15360.d0)
 parameter (J730=100.d0/15360.d0)
 parameter (J731=120.d0/15360.d0)
 parameter (J732=160.d0/15360.d0)
 parameter (J733=320.d0/15360.d0)
 parameter (J740=105.d0/15360.d0)
 parameter (J741=140.d0/15360.d0)
 parameter (J742=280.d0/15360.d0)
 parameter (J750=126.d0/15360.d0)
 parameter (J751=252.d0/15360.d0)
 parameter (J760=231.d0/15360.d0)
 
 parameter (J800=429.d0/34816.d0)
 parameter (J801=462.d0/34816.d0)
 parameter (J802=504.d0/34816.d0)
 parameter (J803=560.d0/34816.d0)
 parameter (J804=640.d0/34816.d0)
 parameter (J805=768.d0/34816.d0)
 parameter (J806=1024.d0/34816.d0)
 parameter (J807=2048.d0/34816.d0)
 parameter (J810=231.d0/34816.d0)
 parameter (J811=252.d0/34816.d0)
 parameter (J812=280.d0/34816.d0)
 parameter (J813=320.d0/34816.d0)
 parameter (J814=384.d0/34816.d0)
 parameter (J815=512.d0/34816.d0)
 parameter (J816=1024.d0/34816.d0)
 parameter (J820=189.d0/34816.d0)
 parameter (J821=210.d0/34816.d0)
 parameter (J822=240.d0/34816.d0)
 parameter (J823=288.d0/34816.d0)
 parameter (J824=284.d0/34816.d0)
 parameter (J825=768.d0/34816.d0)
 parameter (J830=175.d0/34816.d0)
 parameter (J831=200.d0/34816.d0)
 parameter (J832=240.d0/34816.d0)
 parameter (J833=320.d0/34816.d0)
 parameter (J834=640.d0/34816.d0)
 parameter (J840=175.d0/34816.d0)
 parameter (J841=210.d0/34816.d0)
 parameter (J842=280.d0/34816.d0)
 parameter (J843=560.d0/34816.d0)
 parameter (J850=189.d0/34816.d0)
 parameter (J851=252.d0/34816.d0)
 parameter (J852=504.d0/34816.d0)
 parameter (J860=231.d0/34816.d0)
 parameter (J861=462.d0/34816.d0)
 parameter (J870=429.d0/34816.d0)
 
 parameter (J900=6435.d0/622592.d0)
 parameter (J901=6864.d0/622592.d0)
 parameter (J902=7392.d0/622592.d0)
 parameter (J903=8064.d0/622592.d0)
 parameter (J904=8960.d0/622592.d0)
 parameter (J905=10240.d0/622592.d0)
 parameter (J906=12288.d0/622592.d0)
 parameter (J907=16384.d0/622592.d0)
 parameter (J908=32768.d0/622592.d0)
 parameter (J910=3432.d0/622592.d0)
 parameter (J911=3696.d0/622592.d0)
 parameter (J912=4032.d0/622592.d0)
 parameter (J913=4480.d0/622592.d0)
 parameter (J914=5120.d0/622592.d0)
 parameter (J915=6144.d0/622592.d0)
 parameter (J916=8192.d0/622592.d0)
 parameter (J917=16384.d0/622592.d0)
 parameter (J920=2772.d0/622592.d0)
 parameter (J921=3024.d0/622592.d0)
 parameter (J922=3360.d0/622592.d0)
 parameter (J923=3840.d0/622592.d0)
 parameter (J924=4608.d0/622592.d0)
 parameter (J925=6144.d0/622592.d0)
 parameter (J926=12288.d0/622592.d0)
 parameter (J930=2520.d0/622592.d0)
 parameter (J931=2800.d0/622592.d0)
 parameter (J932=3200.d0/622592.d0)
 parameter (J933=3840.d0/622592.d0)
 parameter (J934=5120.d0/622592.d0)
 parameter (J935=10240.d0/622592.d0)
 parameter (J940=2450.d0/622592.d0)
 parameter (J941=2800.d0/622592.d0)
 parameter (J942=3360.d0/622592.d0)
 parameter (J943=4480.d0/622592.d0)
 parameter (J944=8960.d0/622592.d0)
 parameter (J950=2520.d0/622592.d0)
 parameter (J951=3024.d0/622592.d0)
 parameter (J952=4032.d0/622592.d0)
 parameter (J953=8064.d0/622592.d0)
 parameter (J960=2772.d0/622592.d0)
 parameter (J961=3696.d0/622592.d0)
 parameter (J962=7392.d0/622592.d0)
 parameter (J970=3432.d0/622592.d0)
 parameter (J971=6864.d0/622592.d0)
 parameter (J980=6435.d0/622592.d0)
 
 parameter (JA00=12155.d0/1376256.d0)
 parameter (JA01=12870.d0/1376256.d0)
 parameter (JA02=13728.d0/1376256.d0)
 parameter (JA03=14784.d0/1376256.d0)
 parameter (JA04=16128.d0/1376256.d0)
 parameter (JA05=17920.d0/1376256.d0)
 parameter (JA06=20480.d0/1376256.d0)
 parameter (JA07=24576.d0/1376256.d0)
 parameter (JA08=32768.d0/1376256.d0)
 parameter (JA09=65536.d0/1376256.d0)
 parameter (JA10=6435.d0/1376256.d0)
 parameter (JA11=6864.d0/1376256.d0)
 parameter (JA12=7392.d0/1376256.d0)
 parameter (JA13=8064.d0/1376256.d0)
 parameter (JA14=8960.d0/1376256.d0)
 parameter (JA15=10240.d0/1376256.d0)
 parameter (JA16=12288.d0/1376256.d0)
 parameter (JA17=16384.d0/1376256.d0)
 parameter (JA18=32768.d0/1376256.d0)
 parameter (JA20=5148.d0/1376256.d0)
 parameter (JA21=5544.d0/1376256.d0)
 parameter (JA22=6048.d0/1376256.d0)
 parameter (JA23=6720.d0/1376256.d0)
 parameter (JA24=7680.d0/1376256.d0)
 parameter (JA25=9216.d0/1376256.d0)
 parameter (JA26=12288.d0/1376256.d0)
 parameter (JA27=24576.d0/1376256.d0)
 parameter (JA30=4620.d0/1376256.d0)
 parameter (JA31=5040.d0/1376256.d0)
 parameter (JA32=5600.d0/1376256.d0)
 parameter (JA33=6400.d0/1376256.d0)
 parameter (JA34=7680.d0/1376256.d0)
 parameter (JA35=10240.d0/1376256.d0)
 parameter (JA36=20480.d0/1376256.d0)
 parameter (JA40=4410.d0/1376256.d0)
 parameter (JA41=4900.d0/1376256.d0)
 parameter (JA42=5600.d0/1376256.d0)
 parameter (JA43=6720.d0/1376256.d0)
 parameter (JA44=8960.d0/1376256.d0)
 parameter (JA45=17920.d0/1376256.d0)
 parameter (JA50=4410.d0/1376256.d0)
 parameter (JA51=5040.d0/1376256.d0)
 parameter (JA52=6048.d0/1376256.d0)
 parameter (JA53=8064.d0/1376256.d0)
 parameter (JA54=16128.d0/1376256.d0)
 parameter (JA60=4620.d0/1376256.d0)
 parameter (JA61=5544.d0/1376256.d0)
 parameter (JA62=7392.d0/1376256.d0)
 parameter (JA63=14784.d0/1376256.d0)
 parameter (JA70=5148.d0/1376256.d0)
 parameter (JA71=6864.d0/1376256.d0)
 parameter (JA72=13728.d0/1376256.d0)
 parameter (JA80=6435.d0/1376256.d0)
 parameter (JA81=12870.d0/1376256.d0)
 parameter (JA90=12155.d0/1376256.d0)
 
 ! write(*,"(a20,1p3e10.2)") "(serj) y,n,m=",y,n,m
 
 J1=J100
 J2=J200+n*J201+m*J210
 !if(y.lt.9.2925066E-09) then
 !    serj=y*(J1+y*J2)
 ! write(*,"(a20,1pe10.2)") "(serj) J2",serj
 !    return
 !endif
 
 J3=J300+n*(J301+n*J302)+m*(J310+n*J311+m*J320)
 !if(y.lt.4.3667383d-06) then
 !    serj=y*(J1+y*(J2+y*J3))
 ! write(*,"(a20,1pe10.2)") "(serj) J3",serj
 !    return
 !endif
 
 J4=J400+n*(J401+n*(J402+n*J403)) &
     +m*(J410+n*(J411+n*J412)+m*(J420+n*J421+m*J430))
 !if(y.le.9.4990006E-05) then
 !    serj=y*(J1+y*(J2+y*(J3+y*J4)))
 ! write(*,"(a20,1pe10.2)") "(serj) J4",serj
 !    return
 !endif
 
 J5=J500+n*(J501+n*(J502+n*(J503+n*J504))) &
     +m*(J510+n*(J511+n*(J512+n*J513)) &
     +m*(J520+n*(J521+n*J522)+m*(J530+n*J531+m*J540)))
 if(y.le.6.0369310d-04) then
     serj=y*(J1+y*(J2+y*(J3+y*(J4+y*J5))))
 ! write(*,"(a20,1pe10.2)") "(serj) J5",serj
     return
 endif
 
 J6=J600+n*(J601+n*(J602+n*(J603+n*(J604+n*J605)))) &
     +m*(J610+n*(J611+n*(J612+n*(J613+n*J614))) &
     +m*(J620+n*(J621+n*(J622+n*J623)) &
     +m*(J630+n*(J631+n*J632)+m*(J640+n*J641+m*J650))))
 if(y.le.2.0727505d-03) then
     serj=y*(J1+y*(J2+y*(J3+y*(J4+y*(J5+y*J6)))))
 ! write(*,"(a20,1pe10.2)") "(serj) J6",serj
     return
 endif
 
 J7=J700+n*(J701+n*(J702+n*(J703+n*(J704+n*(J705+n*J706))))) &
     +m*(J710+n*(J711+n*(J712+n*(J713+n*(J714+n*J715)))) &
     +m*(J720+n*(J721+n*(J722+n*(J723+n*J724))) &
     +m*(J730+n*(J731+n*(J732+n*J733)) &
     +m*(J740+n*(J741+n*J742)+m*(J750+n*J751+m*J760)))))
 if(y.le.5.0047026d-03) then
     serj=y*(J1+y*(J2+y*(J3+y*(J4+y*(J5+y*(J6+y*J7))))))
 ! write(*,"(a20,1pe10.2)") "(serj) J7",serj
    return
 endif
 
 J8=J800+n*(J801+n*(J802+n*(J803+n*(J804+n*(J805+n*(J806+n*J807)))))) &
     +m*(J810+n*(J811+n*(J812+n*(J813+n*(J814+n*(J815+n*J816))))) &
     +m*(J820+n*(J821+n*(J822+n*(J823+n*(J824+n*J825)))) &
     +m*(J830+n*(J831+n*(J832+n*(J833+n*J834))) &
     +m*(J840+n*(J841+n*(J842+n*J843)) &
     +m*(J850+n*(J851+n*J852)+m*(J860+n*J861+m*J870))))))
 if(y.le.9.6961652d-03) then
     serj=y*(J1+y*(J2+y*(J3+y*(J4+y*(J5+y*(J6+y*(J7+y*J8)))))))
 ! write(*,"(a20,1pe10.2)") "(serj) J8",serj
     return
 endif
 
 J9=J900+n*(J901+n*(J902+n*(J903+n*(J904+n*(J905+n*(J906+n*(J907+n*J908))))))) &
     +m*(J910+n*(J911+n*(J912+n*(J913+n*(J914+n*(J915+n*(J916+n*J917)))))) &
     +m*(J920+n*(J921+n*(J922+n*(J923+n*(J924+n*(J925+n*J926))))) &
     +m*(J930+n*(J931+n*(J932+n*(J933+n*(J934+n*J935)))) &
     +m*(J940+n*(J941+n*(J942+n*(J943+n*J944))) &
     +m*(J950+n*(J951+n*(J952+n*J953)) &
     +m*(J960+n*(J961+n*J962)+m*(J970+n*J971+m*J980)))))))
 if(y.le.1.6220210d-02) then
     serj=y*(J1+y*(J2+y*(J3+y*(J4+y*(J5+y*(J6+y*(J7+y*(J8+y*J9))))))))
 ! write(*,"(a20,1pe10.2)") "(serj) J9",serj
     return
 endif
 
 JA=JA00+n*(JA01+n*(JA02+n*(JA03+n*(JA04+n*(JA05+n*(JA06+n*(JA07+n*(JA08+n*JA09)))))))) &
     +m*(JA10+n*(JA11+n*(JA12+n*(JA13+n*(JA14+n*(JA15+n*(JA16+n*(JA17+n*JA18))))))) &
     +m*(JA20+n*(JA21+n*(JA22+n*(JA23+n*(JA24+n*(JA25+n*(JA26+n*JA27)))))) &
     +m*(JA30+n*(JA31+n*(JA32+n*(JA33+n*(JA34+n*(JA35+n*JA36))))) &
     +m*(JA40+n*(JA41+n*(JA42+n*(JA43+n*(JA44+n*JA45)))) &
     +m*(JA50+n*(JA51+n*(JA52+n*(JA53+n*JA54))) &
     +m*(JA60+n*(JA61+n*(JA62+n*JA63)) &
     +m*(JA70+n*(JA71+n*JA72)+m*(JA80+n*JA81+m*JA90))))))))
 !if(y.lt.2.4482875d-02) then
     serj=y*(J1+y*(J2+y*(J3+y*(J4+y*(J5+y*(J6+y*(J7+y*(J8+y*(J9+y*JA)))))))))
 ! write(*,"(a20,1pe10.2)") "(serj) J10",serj
     return
 !endif
 
 end
 !---------------------------------------------------------------------------
 elemental real(dp) function uatan(t,h)
 !
 ! Universal arctangent function
 !
 ! Reference: Fukushima, T, (2012) J. Comp. Appl. Math., 236, 1961-1975
 !   Precise and fast computation of a general incomplete elliptic integral
 !   of third kind by half and double argument transformations
 !
 real(dp),intent(in) :: t,h !!SJT
 real(dp) :: z,y,x !!SJT
 real(dp) a,r,ri,hold,rold,riold
 real(dp) A3,A5,A7,A9,A11,A13,A15,A17,A19,A21,A23,A25
 !data hold/1.d0/, rold/1.d0/,riold/1.d0/ !!SJT: disable for thread safety
 !save hold,rold,riold !!SJT: disable for thread safety
 parameter (A3=1.d0/3.d0)
 parameter (A5=1.d0/5.d0)
 parameter (A7=1.d0/7.d0)
 parameter (A9=1.d0/9.d0)
 parameter (A11=1.d0/11.d0)
 parameter (A13=1.d0/13.d0)
 parameter (A15=1.d0/15.d0)
 parameter (A17=1.d0/17.d0)
 parameter (A19=1.d0/19.d0)
 parameter (A21=1.d0/21.d0)
 parameter (A23=1.d0/23.d0)
 parameter (A25=1.d0/25.d0)
 !!!!!SJT: replace data statment with manual value assignment for thread safety
 hold=1.d0; rold=1.d0; riold=1.d0
 !!!!!SJT: end replace data statment with manual value assignment for thread safety
 !
 ! write(*,*) "(uatan) t,h=",t,h
 !
 z=-h*t*t
 a=abs(z)
 !
 ! write(*,*) "(uatan) z=",z
 if(a.lt.3.3306691d-16) then
     uatan=t
 ! write(*,*) "0"
 elseif(a.lt.2.3560805d-08) then
     uatan=t*(1.d0+z*A3)
 ! write(*,*) "1"
 elseif(a.lt.9.1939631d-06) then
     uatan=t*(1.d0+z*(A3+z*A5))
 ! write(*,*) "2"
 elseif(a.lt.1.7779240d-04) then
     uatan=t*(1.d0+z*(A3+z*(A5+z*A7)))
 ! write(*,*) "3"
 elseif(a.lt.1.0407839d-03) then
     uatan=t*(1.d0+z*(A3+z*(A5+z*(A7+z*A9))))
 ! write(*,*) "4"
 elseif(a.lt.3.3616998d-03) then
     uatan=t*(1.d0+z*(A3+z*(A5+z*(A7+z*(A9+z*A11)))))
 ! write(*,*) "5"
 elseif(a.lt.7.7408014d-03) then
     uatan=t*(1.d0+z*(A3+z*(A5+z*(A7+z*(A9+z*(A11+z*A13))))))
 ! write(*,*) "6"
 elseif(a.lt.1.4437181d-02) then
     uatan=t*(1.d0+z*(A3+z*(A5+z*(A7+z*(A9+z*(A11+z*(A13+z*A15)))))))
 ! write(*,*) "7"
 elseif(a.lt.2.3407312d-02) then
     uatan=t*(1.d0+z*(A3+z*(A5+z*(A7+z*(A9+z*(A11+z*(A13+z*(A15+z*A17))))))))
 ! write(*,*) "8"
 elseif(a.lt.3.4416203d-02) then
     uatan=t*(1.d0+z*(A3+z*(A5+z*(A7+z*(A9+z*(A11+z*(A13+z*(A15+z*(A17+z*A19)))))))))
 ! write(*,*) "9"
 elseif(z.lt.0.d0) then
     if(abs(h-hold).lt.1.d-16) then
         r=rold; ri=riold
     else
         r=sqrt(h); ri=1.d0/r; hold=h; rold=r; riold=ri
     endif
     uatan=atan(r*t)*ri
 ! write(*,*) "atan"
 elseif(a.lt.4.7138547d-02) then
     uatan=t*(1.d0+z*(A3+z*(A5+z*(A7+z*(A9+z*(A11+z*(A13+z*(A15+z*(A17+z*(A19+z*A21))))))))))
 ! write(*,*) "A"
 elseif(a.lt.6.1227405d-02) then
     uatan=t*(1.d0+z*(A3+z*(A5+z*(A7+z*(A9+z*(A11+z*(A13+z*(A15+z*(A17+z*(A19+z*(A21+z*A23)))))))))))
 ! write(*,*) "B"
 elseif(a.lt.7.6353468d-02) then
     uatan=t*(1.d0+z*(A3+z*(A5+z*(A7+z*(A9+z*(A11+z*(A13+z*(A15+z*(A17+z*(A19+z*(A21+z*(A23+z*A25))))))))))))
 ! write(*,*) "C"
 else
     if(abs(h-hold).lt.1.d-16) then
         r=rold; ri=riold
     else
         r=sqrt(-h); ri=1.d0/r; hold=h; rold=r; riold=ri
     endif
     y=r*t
     x=log((1.d0+y)/(1.d0-y))*0.5d0
     uatan=x*ri
 ! write(*,*) "hyper"
 endif
 !
 ! write(*,*) "(uatan) uatan=",uatan
 !
 return
 end
 
end module xelbdj2_all_routines
