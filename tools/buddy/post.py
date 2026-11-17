import numpy as np, struct, sys, base64, json
from scipy.spatial import cKDTree
fn=sys.argv[1]; d=open(fn,'rb').read(); nv,nt=struct.unpack_from('<II',d,0); o=8
V=np.frombuffer(d,np.float32,nv*3,o).reshape(-1,3); o+=nv*12; T=np.frombuffer(d,np.int32,nt*3,o).reshape(-1,3).copy(); o+=nt*12; C=np.frombuffer(d,np.uint8,nv,o)
# face classes from nearest original faces (several samples per face)
V0=np.load('V.npy'); m=np.load('mesh.npz'); T0=m['CV'].reshape(-1,3); lab=np.load('lab5.npy')
fl0=np.array([np.bincount(r,minlength=5).argmax() for r in lab[T0]]) if False else None
L=lab[T0]; fl0=np.where(L[:,0]==L[:,1],L[:,0],np.where(L[:,1]==L[:,2],L[:,1],np.where(L[:,0]==L[:,2],L[:,0],L[:,0])))
tree=cKDTree(V0[T0].mean(1))
cen=V[T].mean(1); votes=np.zeros((len(T),5))
for w in ([1/3,1/3,1/3],[.6,.2,.2],[.2,.6,.2],[.2,.2,.6]):
    pt=(V[T]*np.array(w)[None,:,None]).sum(1); _,ii=tree.query(pt); np.add.at(votes,(np.arange(len(T)),fl0[ii]),1)
FC=votes.argmax(1)
# smooth normals (area weighted)
fn_=np.cross(V[T[:,1]]-V[T[:,0]],V[T[:,2]]-V[T[:,0]]); N=np.zeros_like(V); 
for j in range(3): np.add.at(N,T[:,j],fn_)
N/=np.linalg.norm(N,axis=1,keepdims=True)+1e-12
# split vertices by face class
key={}; nV=[]; nN=[]; nC=[]; nT=np.zeros_like(T)
for f in range(len(T)):
    for j in range(3):
        k=(T[f,j],FC[f])
        if k not in key: key[k]=len(nV); nV.append(T[f,j]); nC.append(FC[f])
        nT[f,j]=key[k]
src=np.array(nV); P=V[src]; NN=N[src]; CC=np.array(nC)
print('verts',len(P),'tris',len(nT))
np.savez('final.npz',P=P,N=NN,C=CC,T=nT)
