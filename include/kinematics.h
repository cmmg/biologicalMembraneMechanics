/*
  Kirchhoff-Love shell kinematics 
  author: Shiva Rudraraju
*/
#if !defined(KINEMATICS_H_)
#define KINEMATICS_H_

template <class T>
struct  KinematicsStruct{
  //gradient and Hessian maps
  double dXdR[3][2];
  double dXdR2[3][2][2];
  T dxdR[3][2];
  T dxdR2[3][2][2];
  T dxdR_contra[3][2];
  //kinematic quantities
  double A[2][2]; T a[2][2];
  double A_contra[2][2]; T a_contra[2][2];
  double B[2][2]; T b[2][2];
  double B_contra[2][2]; T b_contra[2][2];
  T Gamma[2][2][2];
  //
  T J;
  T normal[3]; double Normal[3];
  //
  T I1;
  double H0;
  T H, Kappa;
};


#undef  __FUNCT__
#define __FUNCT__ "getKinematics"
template <class T>
PetscErrorCode getKinematics(KinematicsStruct<T>& k){
  PetscFunctionBegin;

  //metric tensors A, a
  for(unsigned int i=0; i<2; i++){
    for(unsigned int j=0; j<2; j++){
      k.A[i][j]=0.0;
      k.a[i][j]=0.0;
      for(unsigned int d=0; d<3; d++){
	k.A[i][j]+=k.dXdR[d][i]*k.dXdR[d][j];
	k.a[i][j]+=k.dxdR[d][i]*k.dxdR[d][j];
      }
    }
  }

  //Jacobians
  T J, J_a; double J_A;
  J_A=std::sqrt(k.A[0][0]*k.A[1][1]-k.A[0][1]*k.A[1][0]);
  J_a=std::sqrt(k.a[0][0]*k.a[1][1]-k.a[0][1]*k.a[1][0]);
  if (J_A<=0.0) {std::cout << "negative jacobian\n";  exit(-1);}
  J=J_a/J_A; k.J=J;
  
  //surface normals
  k.normal[0]=(k.dxdR[1][0]*k.dxdR[2][1]-k.dxdR[2][0]*k.dxdR[1][1])/J_a;
  k.normal[1]=(k.dxdR[2][0]*k.dxdR[0][1]-k.dxdR[0][0]*k.dxdR[2][1])/J_a;
  k.normal[2]=(k.dxdR[0][0]*k.dxdR[1][1]-k.dxdR[1][0]*k.dxdR[0][1])/J_a;
  k.Normal[0]=(k.dXdR[1][0]*k.dXdR[2][1]-k.dXdR[2][0]*k.dXdR[1][1])/J_A;
  k.Normal[1]=(k.dXdR[2][0]*k.dXdR[0][1]-k.dXdR[0][0]*k.dXdR[2][1])/J_A;
  k.Normal[2]=(k.dXdR[0][0]*k.dXdR[1][1]-k.dXdR[1][0]*k.dXdR[0][1])/J_A;

  //curvature tensors B, b
  for(unsigned int i=0; i<2; i++){
    for(unsigned int j=0; j<2; j++){
      k.b[i][j]=0.0; k.B[i][j]=0.0;
      for(unsigned int d=0; d<3; d++){
	k.b[i][j]+=k.normal[d]*k.dxdR2[d][i][j];
	k.B[i][j]+=k.Normal[d]*k.dXdR2[d][i][j];
      }
    }
  }

  //determinants of metric tensors and curvature tensor
  T det_a=k.a[0][0]*k.a[1][1]-k.a[0][1]*k.a[1][0];
  T det_b=k.b[0][0]*k.b[1][1]-k.b[0][1]*k.b[1][0];
  double det_A=k.A[0][0]*k.A[1][1]-k.A[0][1]*k.A[1][0];
  
  //contra-variant metric tensors, a_contra, A_contra.
  k.a_contra[0][0]=k.a[1][1]/det_a; k.a_contra[1][1]=k.a[0][0]/det_a;
  k.a_contra[0][1]=-k.a[0][1]/det_a; k.a_contra[1][0]=-k.a[1][0]/det_a;
  k.A_contra[0][0]=k.A[1][1]/det_A; k.A_contra[1][1]=k.A[0][0]/det_A;
  k.A_contra[0][1]=-k.A[0][1]/det_A; k.A_contra[1][0]=-k.A[1][0]/det_A;
  for(unsigned int d=0; d<3; d++){
    k.dxdR_contra[d][0]=k.a_contra[0][0]*k.dxdR[d][0]+k.a_contra[0][1]*k.dxdR[d][1];
    k.dxdR_contra[d][1]=k.a_contra[1][0]*k.dxdR[d][0]+k.a_contra[1][1]*k.dxdR[d][1];
  }

  //contra-variant curvature tensors B_contra, b_contra.
  for(unsigned int i=0; i<2; i++){
    for(unsigned int j=0; j<2; j++){
      k.B_contra[i][j]=0.0;
      k.b_contra[i][j]=0.0;
      for(unsigned int z=0; z<2; z++){
	for(unsigned int l=0; l<2; l++){
	  k.B_contra[i][j]+=k.A_contra[i][z]*k.B[z][l]*k.A_contra[l][j]; //*
	  k.b_contra[i][j]+=k.a_contra[i][z]*k.b[z][l]*k.a_contra[l][j]; //*
	}
      }
    }
  }
  
  //Christoffel symbols
  for(unsigned int i=0; i<2; i++){
    for(unsigned int j=0; j<2; j++){
      for(unsigned int z=0; z<2; z++){
	k.Gamma[i][j][z]=0.0;
	for(unsigned int d=0; d<3; d++){
	  k.Gamma[i][j][z]+=k.dxdR_contra[d][i]*k.dxdR2[d][j][z];
	}
      }
    }
  }

  //invariants of the Green-Lagrange strain, C
  T I1=0.0;
  for(unsigned int i=0; i<2; i++){
    for(unsigned int j=0; j<2; j++){
      I1+=k.A_contra[i][j]*k.a[i][j];
    }
  }
  k.I1=I1;
  
  //mean curvature: H
  T H=0; PetscReal H0=0.0;
  for (unsigned int i=0; i<2; i++){
    for (unsigned int j=0; j<2; j++){
      H+=0.5*k.a[i][j]*k.b_contra[i][j];  //current curvature
      H0+=0.5*k.A[i][j]*k.B_contra[i][j]; //reference curvature
    }
  }
  k.H=H; k.H0=H0;
  
  //Gaussian curvature: Kappa
  k.Kappa=det_b/det_a;
  
  PetscFunctionReturn(0);
}

#endif
