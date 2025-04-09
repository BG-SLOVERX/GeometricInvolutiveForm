read "/mnt/c/Users/Admin/Dropbox/ApproximateCompletion2022/SiyuansVersionsEdits/HybridInvolutiveFormLHPDE.mpl":
read "/mnt/c/Users/Admin/Dropbox/ApproximateCompletion2022/SiyuansVersionsEdits/NumericDiffGeometryTools.mpl":
read "/mnt/c/Users/Admin/Dropbox/ApproximateCompletion2022/SiyuansVersionsEdits/ ApproximateCijk/JoesKam_2_1.dat":

with(RandomTools):
Digits := 15:

targetde := diff(y(x),x,x) + 0.708203932*y(x)*diff(y(x),x) + 5*y(x)^3 + exp(-10((x-1)^2+(y(x)-1)^2));

DimensionInfo := proc(Eqs, randpt)
                local NonlinearDE, DetEqs, rifDetEqs, dimExact, dimExactInfty, dimApprox,GIFinfo, ProlongationOrder, ProlongationTimes;

                NonlinearDE := Eqs;
                DetEqs := DeterminingPDE(NonlinearDE, [xi(x, y), eta(x, y)], integrabilityconditions = false);
                rifDetEqs := rifsimp(DetEqs):
                dimExact := nops(initialdata(rifDetEqs)[Finite]):

                if dimExact = 0 then
                   dimExactInfty := nops(initialdata(rifDetEqs)[Infinite]):
                end if;
    
    
    
    
    # randpt := [rand(1.0..2.0), rand(1.0..2.0)]:
    # dimApprox := GeometricInvolutiveForm(convert(DetEqs, list), [[x, y], [xi, eta]], 10^(-9), 9, randpt):


                # randpt := [1.314, 1.567]:
                GIFinfo := GeometricInvolutiveForm(convert(DetEqs, list), [[x, y], [xi, eta]], 10^(-11), 9, randpt):
                dimApprox := GIFinfo[2]:
                if dimApprox=-1 then
                ProlongationTimes := "skip for now":
                ProlongationOrder := "skip for now":
                else
                ProlongationTimes := GIFinfo[3][1][1]:
                ProlongationOrder := GIFinfo[3][1][2]-1:
                end if;

                #  if dimExact <> dimApprox then
                #     dimExact;
                #     dimApprox;
                #  else
                #     dimApprox;
                #  end if;
                return [dimExact,dimApprox,dimExactInfty, ProlongationTimes, ProlongationOrder, randpt];
                end proc:

# DimensionInfo(de[-2], [1.314,1.567]);

# for j from 1 to 300 do

# fnp := [rand(0.1 .. 5.0)(), rand(0.1 .. 5.0)()]:

# # fnp := [1.0,1.0];

# # for i from -2 to 50 do

# #     # NonlinearDE := de[i];
# #     # # if i = 13 or i = 15 then
# #     # # i;
# #     # # else
# #     # DetEqs := DeterminingPDE(NonlinearDE, [xi(x, y), eta(x, y)], integrabilityconditions = false);
# #     # rifDetEqs := rifsimp(DetEqs):
# #     # dimExact := nops(initialdata(rifDetEqs)[Finite]):
    
    
    
    
# #     # # randpt := [rand(1.0..2.0), rand(1.0..2.0)]:
# #     # # dimApprox := GeometricInvolutiveForm(convert(DetEqs, list), [[x, y], [xi, eta]], 10^(-9), 9, randpt):


# #     # randpt := [1.314, 1.567]:
# #     # dimApprox := GeometricInvolutiveForm(convert(DetEqs, list), [[x, y], [xi, eta]], 10^(-9), 9, randpt)[2]:

# #     # if dimExact <> dimApprox then
# #     #     dimExact;
# #     #     dimApprox;
# #     # else
# #     #     dimApprox;
# #     # end if;
# #     # # end if;
    
# #     DEdimInfo[j][i] := DimensionInfo(de[i], fnp);
# # end do;

# DEdimInfo[j] := DimensionInfo(targetde, fnp):

# end do;
# testtable := [seq(DEdimInfo[i], i=-2..50)]:

DEdiminfoGrid := proc(Eqs, PointSpec)
	local a, b, c, d, nx, nt, dx, dt, i, j, k, DEdimInfoout, DEdimInfotemp;
	userinfo(0, 'NumericDiffGeometryTools', "sys is a system of LHPDE (later extend to poly sys etc)");
	userinfo(0, 'NumericDiffGeometryTools', "vars:  [[indeps], [deps]] 2 indep vars and the corresponding infinitesimal names");
	userinfo(0, 'NumericDiffGeometryTools', "PointSpec is a rectange [[a, b, c, d],[nx,nt]]");
	userinfo(0, 'NumericDiffGeometryTools', "Dmax", "is the max number of prolongations to be used by GIF");

	# if nops(vars[1])<>2 then error `Currently GeometricInvolutiveFormGrid only allows 2 independent variables` end if;

	# if nops(vars[2])<>3 then error `Currently GeometricInvolutiveFormGrid only allows 3 dependent variables` end if;

	# # lprint(`PointSpec=`, PointSpec);

	# if not(type(PointSpec, list)) then error `currently PointSpec must be a list` end if;

	# if nops(PointSpec[1])<>4 or nops(PointSpec[2])<>2 then error `currently must have: nops(PointSpec[1])=4, nops(PointSpec[2])=2` end if;

	# x   := vars[1][1];
	# t   := vars[1][2];
	# xi  := vars[2][1];
	# tau := vars[2][2];
	# phi := vars[2][3];
	
	# a <= x <= b
	# c <= t <= d
	a   := PointSpec[1][1];
	b   := PointSpec[1][2];
	c   := PointSpec[1][3];
	d   := PointSpec[1][4];

	nx  := PointSpec[2][1];
	nt  := PointSpec[2][2];
	dx := (b - a)/nx:
	dt := (d - c)/nt:

    DEdimInfoout := [];

    DEdimInfotemp := [];

	# MList := [ ];

	# TOL := tol;

	# if not(type(TOL,list)) then TOL := [TOL]; end if;
	# TOL := evalf(TOL);

	# if map(type, tol, float) <> [seq(true, i = 1 .. nops(tol))] 
	# 	then error(`You must enter a list of numerical tolerances with`, nops(vars[1]), `entries`);
	# end if;

	# for k to nops(TOL) do
	# 	M  := Matrix(1 .. nx+1, 1 .. nt+1);
		for i from 1 to nx+1 do		
			for j from 1 to nt+1 do
				## lprint("About to enter GIF with i = ", i, "j = ", j):
				# GIF :=     GeometricInvolutiveForm(sys, [[x, t], [xi, tau, phi]], tol, 3, [a+(i-1)*dx, c+(j-1)*dt]);
				#if GIF[2] = -1 then  # i.e. failed to get involutive form
				DEdimInfotemp:= DimensionInfo(Eqs, [a+(i-1/2)*dx, c+(j-1/2)*dt]):
				#end if;
				DEdimInfoout := [op(DEdimInfoout), DEdimInfotemp];
				## lprint("M[",i,j,"]=", M[i,j]):
			end do:
		end do:
	# 	MList := [op(MList), M ];
	# end do:
	return DEdimInfoout;
end proc:


pointspec:= [[0.0, 6.0, 0.0, 6.0], [20,20]]:

DEdimInfo := DEdiminfoGrid(targetde, pointspec);

save DEdimInfo, "targettestfg1gridtol11-0-6-e10.m":