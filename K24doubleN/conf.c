#include "finite.h" 
bool conf(double **pos){
        int l;
        int nb, nd;
        int i, j, k, n;
        double dis, x, y, z, xi, yi, zi;
        int rx = lrint(Nx *0.5);
        int ry = lrint(Ny *0.5);
        int rz = lrint(Nz *0.5);
        double dx = Lx/(Nx-1);
        double dy = Ly/(Ny-1);
        double dz = Lz/(Nz-1);
        double disxy;
        double sinthe, costhe, sinphi, cosphi, omega;
        double Qini[6] = {0};
        double dir[3] = {0};
	int flag = 0;
	double rr;

	        //uniform initial configuration
        if(seed == 0){
                if(!norm_v(init_dir))   return false;
                for(n = 0; n < 6; n ++){
                        Qini[n] = dir2ten(init_dir, n, S);
                }
                for(nd = 0; nd < droplet; nd ++){
                        for(n = 0; n < 6; n ++){
                                Qold[nd * 6 + n] = Qini[n];
                        }
                }
        }

        else if(seed == 1){
                srand(rand_seed);
                for(nd = 0; nd < droplet; nd ++){
                        for(n = 0; n < 3; n ++){
                                dir[n] = (double)rand() / (double)RAND_MAX * 2 - 1;
                        }
                        if(!norm_v(dir))        return false;
                        for(n = 0; n < 6; n ++){
                                Qold[nd * 6 + n] = dir2ten(dir, n, 0.5);
                        }
                }
        }
	//random near particle
	/*
	else if(seed == 10){
		l = 0;
		nd = 0;
                srand(rand_seed);
		if(!norm_v(init_dir))   return false;
		for(k = 0; k < Nz; k++){
			for (j = 0; j < Ny; j++){
				for (i = 0; i < Nx; i++){
					if(drop[l] || boundary[l] || nboundary[l]){
						dir[0] = init_dir[0];			
						dir[1] = init_dir[1];			
						dir[2] = init_dir[2];			
						for(n = 0; n < Np; n++){
							x = i - pos[n][0];
							y = j - pos[n][1];
							z = k - pos[n][2];
							if(x > 0.5 * Nx){
								x -= Nx;
							}
							else if(x < -0.5 * Nx){
								x += Nx;
							}
							if(y > 0.5 * Ny){
								y -= Ny;
							}
							else if(y < -0.5 * Ny){
								y += Ny;
							}
							dis = sqrt(x*x+y*y+z*z);
							if(dis < 2 * Rp){
								for(n = 0; n < 3; n ++){
									dir[n] = (double)rand() / (double)RAND_MAX * 2 - 1;
								}
							}

						}
						if(!norm_v(dir))   return false;
						for(n = 0; n < 6; n ++){
							Qold[nd * 6 + n] = dir2ten(dir, n, S);
						}
						nd ++;
					}
					l ++;
				}
			}
		}
        }
	*/
        //DSS or RSS initial configuration
        else if(seed == 2 || seed == 3){
                l = 0;
		nd = 0;
                for(k = 0; k < Nz; k++){
                        for (j = 0; j < Ny; j++){
                                for (i = 0; i < Nx; i++){
                                        if(drop[l] || boundary[l]||nboundary[l]){
                                                x = (i-rx)*dx;
                                                y = (j-ry)*dy;
                                                z = (k-rz)*dz;
                                                dis = sqrt(x*x+y*y+z*z);
                                                if(seed == 2)   omega = dis * qch;
                                                else
                                                        omega = atan2(y, x) + dis * qch;
                                                disxy = sqrt(x * x + y * y);
                                                if(disxy == 0){
                                                        dir[2] = 1;
                                                        dir[0] = dir[1] = 0;
                                                }
                                                else{
                                                        costhe = z / dis;
                                                        sinthe = disxy / dis;
                                                        cosphi = x / disxy;
                                                        sinphi = y / disxy;
                                                        dir[0] = cos(omega) * costhe * cosphi - sin(omega) * sinphi;
                                                        dir[1] = cos(omega) * costhe * sinphi + sin(omega) * cosphi;
                                                        dir[2] = - cos(omega) * sinthe;
                                                }
                                                if(!norm_v(dir)){
                                                        return false;
                                                }
                                                for (n = 0; n < 6; n++) {
                                                        Qold[nd * 6 + n] = dir2ten(dir, n, S);
                                                }
						nd ++;
                                        }
                                        l ++;
                                }
                        }
                }
        }

        //initial configuration read from Qtensor_ini.out file
        else if(seed == -1){
                if(!norm_v(init_dir))   return false;
                for(n = 0; n < 6; n ++){
                        Qini[n] = dir2ten(init_dir, n, 0.5);
                }
                double a[6] = {0};
                FILE* qtensor;
                qtensor = fopen("Qtensor.bin","rb");
		FILE* grid;
		grid = fopen("grid.bin", "rb");
		int signal;
		if(qtensor == (FILE*)NULL){
			printf("File Qtensor.bin not found.\n");
			return false;
		}
		if(grid == (FILE*)NULL){
			printf("File grid.bin not found.\n");
			return false;
		}
		nd = 0;
		for(l = 0; l < tot; l++){
			fread(&signal, sizeof(int), 1, grid);
			if(signal == 0 || signal == 1){
				fread(a, sizeof(double), 6, qtensor);
				a[5] = - a[0] - a[3];
			}
			else{
				for (n = 0; n < 6; n++) {
					a[n] = Qini[n]; 
				}
			}
			if(drop[l] || boundary[l] || nboundary[l]){
				for (n = 0; n < 6; n++) {
					Qold[nd * 6 + n] = a[n];
				}
				nd ++;
			}
			else{
				for(n = 0; n < 6; n ++){
					Qold[nd * 6 + n] = Qini[n];
				}
			}
		}
		fclose(qtensor);
		fclose(grid);
        }

        //BPI (seed = 4) and BPII (seed = 5)
        else if(seed == 4 || seed == 5){
                double A = 0.2;
                double cst;
                double isq2 = 1.0 / sqrt(2);
                double sq2 = sqrt(2);
                if(seed == 4){
                        cst = 2 * qch * 0.68;
                }
                else{
                        cst = 2 * qch * 0.86;
                }
                l = 0;
		nd = 0;
                for(k = 0; k < Nz; k++){
                        for (j = 0; j < Ny; j++){
                                for (i = 0; i < Nx; i++){
                                        if(drop[l] || boundary[l] || nboundary[l]){
                                                if(seed == 4){
                                                        x = (i - Nx * 0.5) * cst * isq2;
                                                        y = (j - Ny * 0.5) * cst * isq2;
                                                        z = (k - Nz * 0.5) * cst * isq2;
                                                        Qold[nd * 6 + 0] = A * (- sin(y) * cos(x) - sin(x) * cos(z) + 2 * sin(z) * cos(y));
                                                        Qold[nd * 6 + 3] = A * (- sin(z) * cos(y) - sin(y) * cos(x) + 2 * sin(x) * cos(z));
                                                        Qold[nd * 6 + 5] = A * (- sin(x) * cos(z) - sin(z) * cos(y) + 2 * sin(y) * cos(x));
                                                        Qold[nd * 6 + 1] = A * (- sq2 * sin(x) * sin(z) - sq2 * cos(y) * cos(z) + sin(x) * cos(y));
                                                        Qold[nd * 6 + 2] = A * (- sq2 * sin(z) * sin(y) - sq2 * cos(x) * cos(y) + sin(z) * cos(x));
                                                        Qold[nd * 6 + 4] = A * (- sq2 * sin(y) * sin(x) - sq2 * cos(z) * cos(x) + sin(y) * cos(z));
                                                }
                                                else if(seed == 5){
                                                        x = i - Nx * 0.5;
                                                        y = j - Ny * 0.5;
                                                        z = k - Nz * 0.5;
                                                        Qold[nd * 6 + 0] = A * (cos(cst * z) - cos(cst * y));
                                                        Qold[nd * 6 + 1] = A * sin(cst * z);
                                                        Qold[nd * 6 + 2] = A * sin(cst * y);
                                                        Qold[nd * 6 + 3] = A * (cos(cst * x) - cos(cst * z));
                                                        Qold[nd * 6 + 4] = A * sin(cst * x);
                                                        Qold[nd * 6 + 5] = A * (cos(cst * y) - cos(cst * x));
                                                }
						nd ++;
                                        }
                                        l ++;
                                }
                        }
                }
        }
	//seed=6; rotated BPI; 7: rotated BPII
        else if(seed == 6 || seed == 7){
                double A = 0.2;
                double cst;
                double theta = 45 / 180.0 * M_PI;
                double isq2 = 1.0 / sqrt(2);
                double sq2 = sqrt(2);
                if(seed == 6){
                        cst = 2 * qch * 0.68;
                }
                else{
                        cst = 2 * qch * 0.86;
                }
                l = 0;
		nd = 0;
                for(k = 0; k < Nz; k++){
                        for (j = 0; j < Ny; j++){
                                for (i = 0; i < Nx; i++){
                                        if(drop[l] || boundary[l] || nboundary[l]){
                                                if(seed == 6){
                                                        xi = (i - Nx * 0.5) * cst * isq2;
                                                        yi = (j - Ny * 0.5) * cst * isq2;
                                                        zi = (k - Nz * 0.5) * cst * isq2;
							x = xi;
							y = cos(theta) * yi + sin(theta) * zi;
							z = -sin(theta) * yi + cos(theta) * zi;
                                                        Qold[nd * 6 + 0] = A * (- sin(y) * cos(x) - sin(x) * cos(z) + 2 * sin(z) * cos(y));
                                                        Qold[nd * 6 + 3] = A * (- sin(z) * cos(y) - sin(y) * cos(x) + 2 * sin(x) * cos(z));
                                                        Qold[nd * 6 + 5] = A * (- sin(x) * cos(z) - sin(z) * cos(y) + 2 * sin(y) * cos(x));
                                                        Qold[nd * 6 + 1] = A * (- sq2 * sin(x) * sin(z) - sq2 * cos(y) * cos(z) + sin(x) * cos(y));
                                                        Qold[nd * 6 + 2] = A * (- sq2 * sin(z) * sin(y) - sq2 * cos(x) * cos(y) + sin(z) * cos(x));
                                                        Qold[nd * 6 + 4] = A * (- sq2 * sin(y) * sin(x) - sq2 * cos(z) * cos(x) + sin(y) * cos(z));
                                                }
                                                else if(seed == 7){
                                                        xi = i - Nx * 0.5;
                                                        yi = j - Ny * 0.5;
                                                        zi = k - Nz * 0.5;
							x = xi;
							y = cos(theta) * yi + sin(theta) * zi;
							z = -sin(theta) * yi + cos(theta) * zi;
                                                        Qold[nd * 6 + 0] = A * (cos(cst * z) - cos(cst * y));
                                                        Qold[nd * 6 + 1] = A * sin(cst * z);
                                                        Qold[nd * 6 + 2] = A * sin(cst * y);
                                                        Qold[nd * 6 + 3] = A * (cos(cst * x) - cos(cst * z));
                                                        Qold[nd * 6 + 4] = A * sin(cst * x);
                                                        Qold[nd * 6 + 5] = A * (cos(cst * y) - cos(cst * x));
                                                }
						nd ++;
                                        }
                                        l ++;
                                }
                        }
                }
        }
	//helical along z
        else if(seed == 8){
                l = 0;
		nd = 0;
                for(k = 0; k < Nz; k++){
                        for (j = 0; j < Ny; j++){
                                for (i = 0; i < Nx; i++){
                                        if(drop[l] || boundary[l] || nboundary[l]){
						dir[0] = cos(qch * (k - rz));
						dir[1] = sin(qch * (k - rz));
						for(n = 0; n < 6; n ++){
							Qold[nd * 6 + n] = dir2ten(dir, n, S);
						}
						nd ++;
                                        }
                                        l ++;
                                }
                        }
                }
        }
        //deformed DSS or RSS initial configuration in ellip
        else if(seed == 9 || seed == 10){
                l = 0;
		nd = 0;
                for(k = 0; k < Nz; k++){
                        for (j = 0; j < Ny; j++){
                                for (i = 0; i < Nx; i++){
                                        if(drop[l] || boundary[l]||nboundary[l]){
						rr = 2;
                                        //        x = pow(rr, -1/3.0) * (i-rx)*dx;
                                         //       y = pow(rr, -1/3.0) * (j-ry)*dy;
                                          //      z = pow(rr, 2/3.0) * (k-rz)*dz;
                                           //     x = 0.5 * (i-rx)*dx;
                                            //    y = (j-ry)*dy;
                                            //    z = 2 * (k-rz)*dz;
                                                x = 0.5 * pow(rr, -2/3.0) * (i-rx)*dx;
                                                y = pow(rr, 1/3.0) * (j-ry)*dy;
                                                z = 2 * pow(rr, 1/3.0) * (k-rz)*dz;
                                                dis = sqrt(x*x+y*y+z*z);
                                                if(seed == 2)   omega = dis * qch;
                                                else
                                                        omega = atan2(y, x) + dis * qch;
                                                disxy = sqrt(x * x + y * y);
                                                if(disxy == 0){
                                                        dir[2] = 1;
                                                        dir[0] = dir[1] = 0;
                                                }
                                                else{
                                                        costhe = z / dis;
                                                        sinthe = disxy / dis;
                                                        cosphi = x / disxy;
                                                        sinphi = y / disxy;
                                                        dir[0] = cos(omega) * costhe * cosphi - sin(omega) * sinphi;
                                                        dir[1] = cos(omega) * costhe * sinphi + sin(omega) * cosphi;
                                                        dir[2] = - cos(omega) * sinthe;
                                                }
                                                if(!norm_v(dir)){
                                                        return false;
                                                }
                                                for (n = 0; n < 6; n++) {
                                                        Qold[nd * 6 + n] = dir2ten(dir, n, S);
                                                }
						nd ++;
                                        }
                                        l ++;
                                }
                        }
                }
        }
	//dipolar particle
	else if(seed == 11){
		l = 0;
		nd = 0;
		double dip_ini = 2.1;
		if(!norm_v(init_dir))   return false;
		for(k = 0; k < Nz; k++){
			for (j = 0; j < Ny; j++){
				for (i = 0; i < Nx; i++){
					if(drop[l] || boundary[l] || nboundary[l]){
						dir[0] = init_dir[0];			
						dir[1] = init_dir[1];			
						dir[2] = init_dir[2];			
						for(n = 0; n < Np; n++){
							x = i - pos[n][0];
							y = j - pos[n][1];
							z = k - pos[n][2];
							if(x > 0.5 * Nx){
								x -= Nx;
							}
							else if(x < -0.5 * Nx){
								x += Nx;
							}
							if(y > 0.5 * Ny){
								y -= Ny;
							}
							else if(y < -0.5 * Ny){
								y += Ny;
							}
							dis = sqrt(x*x+y*y+z*z);
							dir[0] += dip_ini * Rp * Rp / dis/ dis / dis * x;			
							dir[1] += dip_ini * Rp * Rp / dis/ dis / dis * y;			
							dir[2] += dip_ini * Rp * Rp / dis/ dis / dis * z;			
						}
						if(!norm_v(dir))   return false;
						for(n = 0; n < 6; n ++){
							Qold[nd * 6 + n] = dir2ten(dir, n, S);
						}
						nd ++;
					}
					l ++;
				}
			}
		}
        }
	//twist bipolar particle
	else if(seed == 12){
		l = 0;
		nd = 0;
		if(!norm_v(init_dir))   return false;
		for(k = 0; k < Nz; k++){
			for (j = 0; j < Ny; j++){
				for (i = 0; i < Nx; i++){
					if(drop[l] || boundary[l] || nboundary[l]){
						dir[0] = init_dir[0];			
						dir[1] = init_dir[1];			
						dir[2] = init_dir[2];			
						for(n = 0; n < Np; n++){
							x = i - pos[n][0];
							y = j - pos[n][1];
							z = k - pos[n][2];
							if(x > 0.5 * Nx){
								x -= Nx;
							}
							else if(x < -0.5 * Nx){
								x += Nx;
							}
							if(y > 0.5 * Ny){
								y -= Ny;
							}
							else if(y < -0.5 * Ny){
								y += Ny;
							}
							disxy = sqrt(x*x+y*y);
							dis = sqrt(x*x+y*y + z*z);
							if(fabs(z) > Rp && disxy < Rp && disxy > 0 && fabs(z) < 2 * Rp){
								if(n == 0){
									dir[0] += - 2 * z / fabs(z) *  Rp * Rp / dis/ dis * y / disxy;			
									dir[1] += 2 * z / fabs(z) * Rp * Rp / dis/ dis * x / disxy;			
								}
								else{
									dir[0] += 2 * z / fabs(z) * Rp * Rp / dis/ dis * y / disxy;			
									dir[1] += -2 * z / fabs(z) * Rp * Rp / dis/ dis * x / disxy;			
								}
							}
						}
						if(!norm_v(dir))   return false;
						for(n = 0; n < 6; n ++){
							Qold[nd * 6 + n] = dir2ten(dir, n, S);
						}
						nd ++;
					}
					l ++;
				}
			}
		}
        }
	//radial droplet
        else if(seed == 11){
                l = 0;
		nd = 0;
                for(k = 0; k < Nz; k++){
                        for (j = 0; j < Ny; j++){
                                for (i = 0; i < Nx; i++){
                                        if(drop[l] || boundary[l]||nboundary[l]){
                                                x = (i-rx)*dx;
                                                y = (j-ry)*dy;
                                                z = (k-rz)*dz;
						dir[0] = x;
						dir[1] = y;
						dir[2] = z;
                                                if(!norm_v(dir)){
                                                        dir[0] = 0;
                                                        dir[1] = 0;
                                                        dir[2] = 1;
                                                }
                                                for (n = 0; n < 6; n++) {
                                                        Qold[nd * 6 + n] = dir2ten(dir, n, S);
                                                }
						nd ++;
                                        }
                                        l ++;
                                }
                        }
                }
        }
	//escaped cylinder
	else if(seed == 21){
		l = 0;
		nd = 0;
		double R = rx - 2;
		for(k = 0; k < Nz; k++){
			for (j = 0; j < Ny; j++){
				for (i = 0; i < Nx; i++){
					if(drop[l] || boundary[l] || nboundary[l]){
						x = i - rx;
						y = j - ry;
						dis = sqrt(x*x+y*y);
						dir[0] = x * dis / R;
						dir[1] = y * dis / R;
						dir[2] = sqrt(R * R - dis * dis) / R;
						if(!norm_v(dir))   return false;
						for(n = 0; n < 6; n ++){
							Qold[nd * 6 + n] = dir2ten(dir, n, S);
						}
						nd ++;
					}
					l ++;
				}
			}
		}
        }
        //DSS or RSS initial configuration
        else if(seed == 31){
                l = 0;
		nd = 0;
		double qch0 = qch / 1.1;	
                for(k = 0; k < Nz; k++){
                        for (j = 0; j < Ny; j++){
                                for (i = 0; i < Nx; i++){
                                        if(drop[l] || boundary[l]||nboundary[l]){
                                                x = (i-rx)*dx;
                                                y = (j-ry)*dy;
                                                z = (k-rz)*dz;
                                                dis = sqrt(x*x+y*y+z*z);
                                                omega = atan2(y, x) + dis * qch0;
                                                disxy = sqrt(x * x + y * y);
                                                if(disxy == 0){
                                                        dir[2] = 1;
                                                        dir[0] = dir[1] = 0;
                                                }
                                                else{
                                                        costhe = z / dis;
                                                        sinthe = disxy / dis;
                                                        cosphi = x / disxy;
                                                        sinphi = y / disxy;
                                                        dir[0] = cos(omega) * costhe * cosphi - sin(omega) * sinphi;
                                                        dir[1] = cos(omega) * costhe * sinphi + sin(omega) * cosphi;
                                                        dir[2] = - cos(omega) * sinthe;
                                                }
                                                if(!norm_v(dir)){
                                                        return false;
                                                }
                                                for (n = 0; n < 6; n++) {
                                                        Qold[nd * 6 + n] = dir2ten(dir, n, S);
                                                }
						nd ++;
                                        }
                                        l ++;
                                }
                        }
                }
        }
        else if(seed == 32){
                l = 0;
		nd = 0;
		double qch0 = qch / 1.2;	
                for(k = 0; k < Nz; k++){
                        for (j = 0; j < Ny; j++){
                                for (i = 0; i < Nx; i++){
                                        if(drop[l] || boundary[l]||nboundary[l]){
                                                x = (i-rx)*dx;
                                                y = (j-ry)*dy;
                                                z = (k-rz)*dz;
                                                dis = sqrt(x*x+y*y+z*z);
                                                omega = atan2(y, x) + dis * qch0;
                                                disxy = sqrt(x * x + y * y);
                                                if(disxy == 0){
                                                        dir[2] = 1;
                                                        dir[0] = dir[1] = 0;
                                                }
                                                else{
                                                        costhe = z / dis;
                                                        sinthe = disxy / dis;
                                                        cosphi = x / disxy;
                                                        sinphi = y / disxy;
                                                        dir[0] = cos(omega) * costhe * cosphi - sin(omega) * sinphi;
                                                        dir[1] = cos(omega) * costhe * sinphi + sin(omega) * cosphi;
                                                        dir[2] = - cos(omega) * sinthe;
                                                }
                                                if(!norm_v(dir)){
                                                        return false;
                                                }
                                                for (n = 0; n < 6; n++) {
                                                        Qold[nd * 6 + n] = dir2ten(dir, n, S);
                                                }
						nd ++;
                                        }
                                        l ++;
                                }
                        }
                }
        }
        else if(seed == 33){
                l = 0;
		nd = 0;
		double qch0 = qch / 1.3;	
                for(k = 0; k < Nz; k++){
                        for (j = 0; j < Ny; j++){
                                for (i = 0; i < Nx; i++){
                                        if(drop[l] || boundary[l]||nboundary[l]){
                                                x = (i-rx)*dx;
                                                y = (j-ry)*dy;
                                                z = (k-rz)*dz;
                                                dis = sqrt(x*x+y*y+z*z);
                                                omega = atan2(y, x) + dis * qch0;
                                                disxy = sqrt(x * x + y * y);
                                                if(disxy == 0){
                                                        dir[2] = 1;
                                                        dir[0] = dir[1] = 0;
                                                }
                                                else{
                                                        costhe = z / dis;
                                                        sinthe = disxy / dis;
                                                        cosphi = x / disxy;
                                                        sinphi = y / disxy;
                                                        dir[0] = cos(omega) * costhe * cosphi - sin(omega) * sinphi;
                                                        dir[1] = cos(omega) * costhe * sinphi + sin(omega) * cosphi;
                                                        dir[2] = - cos(omega) * sinthe;
                                                }
                                                if(!norm_v(dir)){
                                                        return false;
                                                }
                                                for (n = 0; n < 6; n++) {
                                                        Qold[nd * 6 + n] = dir2ten(dir, n, S);
                                                }
						nd ++;
                                        }
                                        l ++;
                                }
                        }
                }
        }
        else if(seed == 34){
                l = 0;
		nd = 0;
		double qch0 = qch / 1.4;	
                for(k = 0; k < Nz; k++){
                        for (j = 0; j < Ny; j++){
                                for (i = 0; i < Nx; i++){
                                        if(drop[l] || boundary[l]||nboundary[l]){
                                                x = (i-rx)*dx;
                                                y = (j-ry)*dy;
                                                z = (k-rz)*dz;
                                                dis = sqrt(x*x+y*y+z*z);
                                                omega = atan2(y, x) + dis * qch0;
                                                disxy = sqrt(x * x + y * y);
                                                if(disxy == 0){
                                                        dir[2] = 1;
                                                        dir[0] = dir[1] = 0;
                                                }
                                                else{
                                                        costhe = z / dis;
                                                        sinthe = disxy / dis;
                                                        cosphi = x / disxy;
                                                        sinphi = y / disxy;
                                                        dir[0] = cos(omega) * costhe * cosphi - sin(omega) * sinphi;
                                                        dir[1] = cos(omega) * costhe * sinphi + sin(omega) * cosphi;
                                                        dir[2] = - cos(omega) * sinthe;
                                                }
                                                if(!norm_v(dir)){
                                                        return false;
                                                }
                                                for (n = 0; n < 6; n++) {
                                                        Qold[nd * 6 + n] = dir2ten(dir, n, S);
                                                }
						nd ++;
                                        }
                                        l ++;
                                }
                        }
                }
        }
	return true;
}
