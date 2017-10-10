#!/usr/bin/awk

 
### 128bit xorshift
function xorshift(){
	t = xor(x, lshift(x,11))
	x = y; y = z; z = w;
	return w = xor(xor(w,rshift(w,19)), xor(t,rshift(t,8)))
}

function read_file(array_FILE,	array) {
    while ((getline < array_FILE) > 0){
		for (id=1; id<=3; id++){
			num=id+1
        	array[$1,id] = $num
		}
    }
    close (array_FILE)
}

function calc_force_bond(array, force){
	for (i=1; i<=N_atom-1; i++) {
		for (id=1; id<=3; id++){
			vec[i,id]=array[i+1,id]-array[i,id]
		}
		Lbond=sqrt(vec[i,1]^2+vec[i,2]^2+vec[i,3]^2)
		dbond=Lbond-s_bond[i]
		ddbond=dbond^2
		f=(k_bond+k_bond*2*ddbond)*(-2*dbond/Lbond)
		for (id=1; id<=3; id++){
			f_bond[i,id]=f*vec[i,id]
			force[i,id,1]-=f_bond[i,id]
			force[i+1,id,1]+=f_bond[i,id]
		}
	}
}

function calc_force_angl(array, force){
	for (i=1; i<=N_atom-2; i++) {
		for (id=1; id<=3; id++){
			vec1[i,id]=array[i+1,id]-array[i,id]
			vec2[i,id]=array[i+2,id]-array[i+1,id]
		}
		vec1d=vec1[i,1]^2+vec1[i,2]^2+vec1[i,3]^2
		vec2d=vec2[i,1]^2+vec2[i,2]^2+vec2[i,3]^2
		naiseki=vec1[i,1]*vec2[i,1]+vec1[i,2]*vec2[i,2]+vec1[i,3]*vec2[i,3]
		cos_thita=-naiseki/sqrt(vec1d*vec2d)
		if (cos_thita > 1) {
			cos_thita=1
		} else if (cos_thita < -1){
			cos_thita=-1
		}
		sin_thita=sqrt(1-cos_thita^2)
		thita=atan2(sin_thita,cos_thita)
		dba=thita-s_angl[i]
		t3=vec1d*vec2d-naiseki^2
		if (t3 < 1) {
			t3=1
		}
		f=2*k_angl*dba/sqrt(t3)
		f-=k_angl/sqrt(vec1d*vec2d)
		for (id=1; id<=3; id++){
			f_angl21=f*(vec1[i,id]*(naiseki/vec1d)-vec2[i,id])
			f_angl32=f*(vec2[i,id]*(naiseki/vec2d)-vec1[i,id])
			force[i,id,1]-=f_angl21
			force[i+1,id,1]+=f_angl21-f_angl32
			force[i+2,id,1]+=f_angl32
		}
	}
}

function calc_force_exv(array, force){
    cutoff2=cutoff^2; d_exv2=d_exv^2
    for (i=1; i<=N_atom-2; i++){
        for (j=i+2; j<=N_atom; j++){
		    for (id=1; id<=3; id++){
                vec[j,id]=array[i,id]-array[j,id]
            }
            dist2=vec[j,1]^2+vec[j,2]^2+vec[j,3]^2
            if (dist2 > cutoff2) {
                break
            }
            f=k_exv*(d_exv2/dist2)^7
		    for (id=1; id<=3; id++){
                f_exv=f*vec[j,id]
                force[i,id,1]+=f_exv
                force[j,id,1]-=f_exv
            }
        }
    }
}


function calc_xyz(array, velo, force, rforce){
	for (i=1; i<=N_atom; i++){
		for (id=1; id<=3; id++){
			array[i,id]+=velo[i,id]*dt*(1-(gamma*dt)/2)
			array[i,id]+=dt^2*0.5*(force[i,id,0]/m+rforce[i,id,0])
		}
	}
}

function calc_velo(velo, force, rforce){
	for (i=1; i<=N_atom; i++){
		for (id=1; id<=3; id++){
			rforce[i,id,1]=rand_force()
			velo[i,id]=velo[i,id]*(1-(gamma*dt)*0.5)*(1-(gamma*dt)*0.5+(gamma*dt)^2*0.25)
			velo[i,id]+=dt*0.5*(1-(gamma*dt)*0.5)*((force[i,id,0]-force[i,id,1])/m+rforce[i,id,0]+rforce[i,id,1])
		}
	}
}

function rand_force(r){   ### 値のオーダーがforceの0.1倍も小さいんだけどいいの？あと、３方向への配分とかは？
	max_score=2^53
	ran1 = xorshift()/max_score
	ran2 = xorshift()/max_score
	r = sqrt(-2.0*log(ran1))*cos(2.0*atan2(0,-0)*ran2)
	r=r*2*gamma*kb*T/(m*dt)
    return r
}

function calc_energy_bond(array){
    e_bond=0
	for (i=1; i<=N_atom-1; i++) {
		Lbond=sqrt((array[i+1,1]-array[i,1])^2+(array[i+1,2]-array[i,2])^2+(array[i+1,3]-array[i,3])^2)
		ddbond=(Lbond-s_bond[i])^2
		e_bond+=(k_bond+k_bond*ddbond)*ddbond
	}
}

function calc_energy_angl(array){
    e_angl=0
	for (i=1; i<=N_atom-2; i++) {
		for (id=1; id<=3; id++){
			vec1[i,id]=array[i+1,id]-array[i,id]
			vec2[i,id]=array[i+2,id]-array[i+1,id]
		}
		vec1d=vec1[i,1]^2+vec1[i,2]^2+vec1[i,3]^2
		vec2d=vec2[i,1]^2+vec2[i,2]^2+vec2[i,3]^2
		naiseki=vec1[i,1]*vec2[i,1]+vec1[i,2]*vec2[i,2]+vec1[i,3]*vec2[i,3]
		cos_thita=-naiseki/sqrt(vec1d*vec2d)
		if (cos_thita > 1) {
			cos_thita=1
		} else if (cos_thita < -1){
			cos_thita=-1
		}
		sin_thita=sqrt(1-cos_thita^2)
		thita=atan2(sin_thita,cos_thita)
		dba=thita-s_angl[i]
		e_angl+=k_angl*dba^2+k_angl*(cos_thita+1)
	}
}

function calc_energy_exv(array){
    cutoff2=cutoff^2; d_exv2=d_exv^2; e_exv=0
    for (i=1; i<=N_atom-2; i++){
        for (j=i+2; j<=N_atom; j++){
		    for (id=1; id<=3; id++){
                vec[j,id]=array[i,id]-array[j,id]
            }
            dist2=vec[j,1]^2+vec[j,2]^2+vec[j,3]^2
            if (dist2 > cutoff2) {
                break
            }
            e_exv+=k_exv*(d_exv2/dist2)^6
        }
    }
}

function output_trajectory(array){
    print "MODEL" >> "trajectory.movie"
    print "<<<<" >> "trajectory.movie"
    for (i=1; i<=N_atom; i++){
        printf("%4s%7s%6s%3s%2s%4s%4s%8.3f%8.3f%8.3f%12s\n","ATOM",i,"  CA  ",HU_arr[i]," A",i,"    ",array[i,1],array[i,2],array[i,3],"  0.00  0.00") >> "trajectory.movie"
    }
    print ">>>>" >> "trajectory.movie"
    print "ENDMDL" >> "trajectory.movie"
}

function output_energy(array){
    calc_energy_bond(array)
    calc_energy_angl(array)
    calc_energy_exv(array)
    total=e_bond+e_angl+e_exv
    print loop,e_bond,e_angl,e_exv,total >> "energy.ts"
}

function shuffle_num(HU_bind_arr){
	for (i=1; i<=N_atom; i++){
		HU_bind_arr[i]=i
	}
	for (i=N_atom; i >= 1; i--) {
		j = int((i+1)*rand())
		if (j == 0) {
			j = 1
		}
		tmp = HU_bind_arr[i]
		HU_bind_arr[i] = HU_bind_arr[j]
		HU_bind_arr[j] = tmp
	}
}

BEGIN{
srand()
x = 123456789; y = 362436069; z = 521288629; w = 88675123; ### for xorshift
w=w*int(rand()*10)   ### 32bit整数をきちんとrand()から出力するには？

dt=0.1; gamma=0.25; m=100; kb=0.002; T=300.0   ### kb*T=0.6
k_bond=10.0; k_angl=1.0; k_exv=0.6
cutoff=2.0; d_exv=1.0
loop=0

read_file(ARGV[1], array)
N_atom=length(array)/3

shuffle_num(HU_bind_arr)
HU_num=1

for (i=1; i<=N_atom; i++){
	for (id=1; id<=3; id++){
		velo[i,id]=0
	}
	s_bond[i]=1.0
	s_angl[i]=180.0*(atan2(0,-0)/180)
	HU_arr[i]="XXX"
}

output_trajectory(array)
output_energy(array)

calc_force_bond(array, force)
calc_force_angl(array, force)
calc_force_exv(array, force)
for (i=1; i<=N_atom; i++){
	for (id=1; id<=3; id++){
        force[i,id,0]=force[i,id,1]
		rforce[i,id,0]=rand_force()
	}
}

### loop
for (loop=1; loop<=1000000; loop++){
	calc_xyz(array, velo, force, rforce)
	calc_force_bond(array, force)
	calc_force_angl(array, force)
    calc_force_exv(array, force)
	calc_velo(velo, force, rforce)

    if (loop%100 == 0) {
        output_trajectory(array)
    }

    if (loop%1000 == 0) {
        output_energy(array)
	}

	if (loop%10000 == 0 && loop <= 300000) {
		HU_bind_num=HU_bind_arr[HU_num]		
		s_angl[HU_bind_num]=100.0*(atan2(0,-0)/180)
		HU_arr[HU_bind_num+1]="AAA"
		HU_num+=1
	}

    for (i=1; i<=N_atom; i++){
        for (id=1; id<=3; id++){
            force[i,id,0]=force[i,id,1]
            rforce[i,id,0]=rforce[i,id,1]
            force[i,id,1]=0.0
            rforce[i,id,1]=0.0
        }
    }
}
### loop

}
