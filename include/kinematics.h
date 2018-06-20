/*
  Kirchhoff-Love shell kinematics 
  author: Shiva Rudraraju
*/
#if !defined(KINEMATICS_H_)
#define KINEMATICS_H_

template <class T>
typedef struct {
  //gradient and Hessian maps
  double dXdR[3][2][2];
  double dXdR2[3][2][2];
  T dxdR[3][2];
  T dxdR2[3][2][2];

  //kinematic quantities
  double A[2][2]; T a[2][2];
  double A_contra[2][2]; T a_contra[2][2];
  double B[2][2]; T b[2][2];
  double B_contra[2][2]; T b_contra[2][2];
  T Gamma[2][2][2]
  //
  T J;
  T normal[3]; double Normal[3];
  //
  T I1;
  T H, Kappa;
} KinematicsStruct;


#undef  __FUNCT__
#define __FUNCT__ "getKinematics"
template <class T>
PetscErrorCode getKinematics(void* _kinematics){
  PetscFunctionBegin;
  KinematicsStruct<T> *k = (KinematicsStruct<T> *) _kinematics;
  
  //metric tensors A, a
  double (*A)[2] = (double (*)[2]) k->A;
  T (*a)[2] = (T (*)[2]) k->a;
  for(unsigned int i=0; i<2; i++){
    for(unsigned int j=0; j<2; j++){
      A[i][j]=0.0;
      a[i][j]=0.0;
      for(unsigned int d=0; d<3; d++){
	A[i][j]+=dXdR[d][i]*dXdR[d][j];
	a[i][j]+=dxdR[d][i]*dxdR[d][j];
      }
    }
  }

  //Jacobians
  T J, J_a; double J_A;
  J_A=std::sqrt(A[0][0]*A[1][1]-A[0][1]*A[1][0]);
  J_a=std::sqrt(a[0][0]*a[1][1]-a[0][1]*a[1][0]);
  if (J_A<=0.0) {std::cout << "negative jacobian\n";  exit(-1);}
  J=J_a/J_A; k->J=J;
  
  //surface normals
  double (*Normal) = (double (*)) k->Normal;
  T (*normal) = (T (*)) k->normal;
  normal[0]=(dxdR[1][0]*dxdR[2][1]-dxdR[2][0]*dxdR[1][1])/J_a;
  normal[1]=(dxdR[2][0]*dxdR[0][1]-dxdR[0][0]*dxdR[2][1])/J_a;
  normal[2]=(dxdR[0][0]*dxdR[1][1]-dxdR[1][0]*dxdR[0][1])/J_a;
  Normal[0]=(dXdR[1][0]*dXdR[2][1]-dXdR[2][0]*dXdR[1][1])/J_A;
  Normal[1]=(dXdR[2][0]*dXdR[0][1]-dXdR[0][0]*dXdR[2][1])/J_A;
  Normal[2]=(dXdR[0][0]*dXdR[1][1]-dXdR[1][0]*dXdR[0][1])/J_A;
  
  //curvature tensors B, b
  double (*B)[2] = (double (*)[2]) k->B;
  T (*b)[2] = (T (*)[2]) k->b;
  for(unsigned int i=0; i<2; i++){
    for(unsigned int j=0; j<2; j++){
      b[i][j]=0.0; B[i][j]=0.0;
      for(unsigned int d=0; d<3; d++){
	b[i][j]+=normal[d]*dxdR2[d][i][j];
	B[i][j]+=Normal[d]*dXdR2[d][i][j];
      }
    }
  }

  //determinants of metric tensors and curvature tensor
  T det_a=a[0][0]*a[1][1]-a[0][1]*a[1][0];
  T det_b=b[0][0]*b[1][1]-b[0][1]*b[1][0];
  double det_A=A[0][0]*A[1][1]-A[0][1]*A[1][0];
  
  //contra-variant metric tensors, a_contra, A_contra.
  double (*A_contra)[2] = (double (*)[2]) k->A_contra;
  T (*a_contra)[2] = (T (*)[2]) k->a_contra;
  T dxdR_contra[3][2];
  a_contra[0][0]=a[1][1]/det_a; a_contra[1][1]=a[0][0]/det_a;
  a_contra[0][1]=-a[0][1]/det_a; a_contra[1][0]=-a[1][0]/det_a;
  A_contra[0][0]=A[1][1]/det_A; A_contra[1][1]=A[0][0]/det_A;
  A_contra[0][1]=-A[0][1]/det_A; A_contra[1][0]=-A[1][0]/det_A;
  for(unsigned int d=0; d<3; d++){
    dxdR_contra[d][0]=a_contra[0][0]*dxdR[d][0]+a_contra[0][1]*dxdR[d][1];
    dxdR_contra[d][1]=a_contra[1][0]*dxdR[d][0]+a_contra[1][1]*dxdR[d][1];
  }
  
  //contra-variant curvature tensors B_contra, b_contra.
  double (*B_contra)[2] = (double (*)[2]) k->B_contra;
  T (*b_contra)[2] = (T (*)[2]) k->b_contra;
  for(unsigned int i=0; i<2; i++){
    for(unsigned int j=0; j<2; j++){
      B_contra[i][j]=0.0;
      b_contra[i][j]=0.0;
      for(unsigned int k=0; k<2; k++){
	for(unsigned int l=0; l<2; l++){
	  B_contra[i][j]+=A_contra[i][k]*B[k][l]*A_contra[l][j]; //*
	  b_contra[i][j]+=a_contra[i][k]*b[k][l]*a_contra[l][j]; //*
	}
      }
    }
  }
  
  //Christoffel symbols
  T (*Gamma)[2][2] = (T (*)[2][2]) k->Gamma;
  for(unsigned int i=0; i<2; i++){
    for(unsigned int j=0; j<2; j++){
      for(unsigned int k=0; k<2; k++){
	Gamma[i][j][k]=0.0;
	for(unsigned int d=0; d<3; d++){
	  Gamma[i][j][k]+=dxdR_contra[d][i]*dxdR2[d][j][k];
	}
      }
    }
  }

  //invariants of the Green-Lagrange strain, C
  T I1=0.0;
  for(unsigned int i=0; i<2; i++){
    for(unsigned int j=0; j<2; j++){
      I1+=A_contra[i][j]*a[i][j];
    }
  }
  k->I1=I1;
  
  //mean curvature: H
  T H=0;
  for (unsigned int i=0; i<2; i++){
    for (unsigned int j=0; j<2; j++){
      H+=0.5*a[i][j]*b_contra[i][j];  //current curvature
      H0+=0.5*A[i][j]*B_contra[i][j]; //reference curvature
    }
  }
  k->H=H;
  
  //Gaussian curvature: Kappa
  T Kappa=det_b/det_a;
  k->Kappa=Kappa;
  
  PetscFunctionReturn(0);
}

#endif
