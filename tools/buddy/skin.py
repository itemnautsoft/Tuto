import numpy as np, base64, json, sys
f=np.load('final.npz'); P=f['P'].copy(); N=f['N']; C=f['C']; T=f['T']
S=1.55
ax=np.abs(P[:,0]); sg=np.sign(P[:,0])
armband=(P[:,1]>0.72*S)&(P[:,1]<0.94*S)
# shorten the arms to Tuto's reach: tube squeezed, mitten moved in rigidly
t0,t1,dl=0.17*S,0.40*S,-0.15
k=(t1+dl-t0)/(t1-t0)
m=armband&(ax>t0)
nx=np.where(ax<t1, t0+(ax-t0)*k, ax+dl)
P[m,0]=(sg*nx)[m]
ss=lambda a,b,x: np.clip((x-a)/(b-a),0,1)**2*(3-2*np.clip((x-a)/(b-a),0,1))
B=['body','head','spin','shL','elL','shR','elR','hipL','knL','hipR','knR']
W=np.zeros((len(P),len(B)))
x,y=P[:,0],P[:,1]; ax=np.abs(x)
wh=ss(1.42,1.56,y); wsp=ss(2.64,2.70,y)
W[:,1]=wh*(1-wsp); W[:,2]=wh*wsp
rest=1-wh
isArm=(y>1.08)&(y<1.46)&(ax>0.22)
a1=ss(0.22,0.33,ax); a2=ss(0.36,0.47,ax)
L=x<0
legW=ss(1.02,0.80,y); kk=ss(0.60,0.46,y); sL=ss(0.04,-0.04,x)
for i in range(len(P)):
    if isArm[i]:
        W[i,0]+=rest[i]*(1-a1[i]); s,e=(3,4) if L[i] else (5,6)
        W[i,s]+=rest[i]*a1[i]*(1-a2[i]); W[i,e]+=rest[i]*a1[i]*a2[i]
    else:
        lw=rest[i]*legW[i]; W[i,0]+=rest[i]-lw
        for side,hb,kb in ((sL[i],7,8),(1-sL[i],9,10)):
            W[i,hb]+=lw*side*(1-kk[i]); W[i,kb]+=lw*side*kk[i]
idx=np.argsort(-W,axis=1)[:,:4]; w=np.take_along_axis(W,idx,1); w/=w.sum(1,keepdims=True)
wq=np.round(w*255).astype(int); wq[:,0]+=255-wq.sum(1)
print('bone usage',[(B[b],int((idx[:,0]==b).sum())) for b in range(len(B))])
Q=np.round(P*10000).astype(np.int16)
b64=lambda a: base64.b64encode(a.tobytes()).decode()
asset=dict(n=len(P),p=b64(Q),nr=b64(np.round(N*127).astype(np.int8)),c=b64(C.astype(np.uint8)),i=b64(T.astype(np.uint16)),si=b64(idx.astype(np.uint8)),sw=b64(wq.astype(np.uint8)))
s='const BUDDY = '+json.dumps(asset)+';\n'
open('buddy_asset.js','w').write(s); print('asset bytes',len(s))
