"""Equivalent implementations of the reference MAC operators.

Fuse axial/angular Fourier products; reuse sparse Stokes patterns; solve
signed axial frequency pairs together. No PDE/model/boundary change.
"""
import numpy as np
from scipy.fft import fft,ifft,fftn,ifftn
from scipy.sparse import diags,bmat
from scipy.sparse.linalg import splu
from evolve_full_mac import MAC


class FastMAC(MAC):
    def nonlinear(self,U,project=True):
        ah=fft(U,axis=-1,norm='forward')
        uc=np.stack([self.gather(a) for a in ah])
        wc=np.empty_like(uc)
        for j,m in enumerate(self.modes):
            u=uc[j]
            wc[j]=np.stack((1j*m*u[2]/self.r[:,None]-1j*self.kz*u[1],
                1j*self.kz*u[0]-self.Dc_even@u[2],
                self.Dc_odd@u[1]+u[1]/self.r[:,None]-1j*m*u[0]/self.r[:,None]))
        nwork=self.nzpad if self.dealias=='padding' else self.nz
        positive_indices=self.ki%nwork;negative_indices=(-self.ki)%nwork
        cell_force=np.zeros_like(uc)
        for i in range(0,self.nr,24):
            sl=slice(i,min(i+24,self.nr));shape=(3,sl.stop-sl.start,nwork,self.nphi)
            us=np.zeros(shape,complex);ws=np.zeros(shape,complex)
            for j in range(self.J+1):
                us[...,j][:,:,positive_indices]=uc[j,:,sl]
                ws[...,j][:,:,positive_indices]=wc[j,:,sl]
                if j:
                    # Reality reverses BOTH Fourier indices before conjugation.
                    us[...,-j][:,:,negative_indices]=uc[j,:,sl].conj()
                    ws[...,-j][:,:,negative_indices]=wc[j,:,sl].conj()
            up=ifftn(us,axes=(-2,-1),norm='forward').real
            wp=ifftn(ws,axes=(-2,-1),norm='forward').real
            cross=np.stack((up[1]*wp[2]-up[2]*wp[1],up[2]*wp[0]-up[0]*wp[2],up[0]*wp[1]-up[1]*wp[0]))
            result=fftn(cross,axes=(-2,-1),norm='forward')
            for j in range(self.J+1):cell_force[j,:,sl]=result[...,j][:,:,positive_indices]
        fh=np.stack([self.scatter(a) for a in cell_force])
        fh[:,:,~self.keep]=0
        if project:
            for j in range(self.J+1):fh[j]=self.project_hat(fh[j],j)
        result=ifft(fh,axis=-1,norm='forward');result[0]=result[0].real
        return result

    def _stokes_template(self,j,alpha,zero):
        if not hasattr(self,'pattern_cache'):self.pattern_cache={}
        key=(j,alpha,zero)
        if key in self.pattern_cache:return self.pattern_cache[key]
        n=self.nr;m=self.modes[j]
        D=bmat([[self.Divr,diags(1j*m/self.r),1j*diags(np.ones(n))]],format='csc')
        BM=-D.conj().T@diags(self.V)
        full=bmat([[diags(self.mass)+alpha*self.K[j],BM],[BM.conj().T,None]],format='csc')
        perm=[]
        for i in range(n):
            perm.extend((n-1+i,2*n-1+i))
            if i<n-1:perm.append(i)
            if not(zero and i==n-1):perm.append(self.nv+i)
        perm=np.array(perm);matrix=full[perm][:,perm].tocsc();matrix.sort_indices()
        inverse={int(original):new for new,original in enumerate(perm)}
        def position(row,col):
            r=inverse[row];c=inverse[col]
            begin,end=matrix.indptr[c:c+2]
            index=begin+np.searchsorted(matrix.indices[begin:end],r)
            assert index<end and matrix.indices[index]==r
            return index
        diagonal=np.array([position(i,i) for i in range(self.nv)])
        linear=[]
        for i in range(n):
            if zero and i==n-1:continue
            zi=2*n-1+i;pi=self.nv+i
            linear.extend((position(zi,pi),position(pi,zi)))
        value=(matrix,perm,diagonal,np.array(linear))
        self.pattern_cache[key]=value
        return value

    def _stokes_factor(self,j,kindex,alpha):
        key=(j,kindex,alpha)
        if key in self.stokes_cache:return self.stokes_cache[key]
        matrix,perm,diagonal,linear=self._stokes_template(j,alpha,j==0 and kindex==0)
        k=2*np.pi*kindex/self.L
        a=matrix.copy();a.data[diagonal]+=alpha*k*k*self.mass
        a.data[linear]*=k
        value=(splu(a,permc_spec='NATURAL'),perm)
        self.stokes_cache[key]=value
        return value

    def diffuse(self,U,h):
        alpha=h*self.nu/2;out=np.zeros_like(U)
        for j in range(self.J+1):
            ah=fft(U[j],axis=-1)
            rhs=self.mass[:,None]*ah-alpha*(self.K[j]@ah+self.mass[:,None]*self.kz**2*ah)
            vh=np.zeros_like(ah)
            for kindex in np.flatnonzero(self.keep&(self.ki>=0)):
                cols=[int(kindex)] if kindex==0 else [int(kindex),self.nz-int(kindex)]
                lu,perm=self._stokes_factor(j,int(kindex),alpha)
                b=np.zeros((self.nv+self.nr,len(cols)),complex);b[:self.nv]=rhs[:,cols]
                if kindex:b[self.slz,1]*=-1
                solution=np.zeros_like(b)
                solution[perm]=lu.solve(np.asfortranarray(b[perm]))
                if kindex:solution[self.slz,1]*=-1
                vh[:,cols]=solution[:self.nv]
            out[j]=ifft(vh,axis=-1)
        out[0]=out[0].real
        return out
