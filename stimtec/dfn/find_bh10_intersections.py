import vtk, numpy as np
from vtk.util.numpy_support import vtk_to_numpy as v2n
from scipy.spatial import cKDTree
np.set_printoptions(precision=4, suppress=True, linewidth=200)
def rd(f):
    r=vtk.vtkXMLUnstructuredGridReader(); r.SetFileName(f); r.Update(); return r.GetOutput()
m=rd('mesh.vtu'); b=rd('BH10.vtu')
P=v2n(m.GetPoints().GetData()); C=v2n(m.GetCells().GetConnectivityArray()).reshape(-1,3)
mat=v2n(m.GetCellData().GetArray('MaterialIDs'))
B=v2n(b.GetPoints().GetData())
# BH10 points vs mesh nodes
d,i=cKDTree(P).query(B); print('BH10 pts nearest mesh node dist:',d.round(3))
bulk=v2n(b.GetPointData().GetArray('bulk_node_ids')); print('bulk_node_ids coords dist:',np.linalg.norm(P[bulk.astype(int)]-B,axis=1).round(2))
# polyline: start, end, direction
A=B[0]; E=B[-1]; L=np.linalg.norm(E-A); u=(E-A)/L
print('BH10 length',L,'dir',u, 'collinear dev', np.abs(np.cross(B-A,u)).max())
# segment-triangle intersection (Moller-Trumbore) for whole line A->E
v0,v1,v2=P[C[:,0]],P[C[:,1]],P[C[:,2]]
e1=v1-v0; e2=v2-v0; h=np.cross(u,e2); a=(e1*h).sum(1)
ok=np.abs(a)>1e-12; f=1/np.where(ok,a,1); s=A-v0
uu=f*(s*h).sum(1); q=np.cross(s,e1); vv=f*(u*q).sum(1); t=f*(e2*q).sum(1)
tol=1e-9
hit=ok&(uu>=-tol)&(vv>=-tol)&(uu+vv<=1+tol)&(t>=-1e-6)&(t<=L+1e-6)
hits=np.where(hit)[0]
print('n hit triangles',len(hits))
res=[]
for c in hits:
    X=A+t[c]*u; res.append((t[c],c,mat[c],X))
res.sort(key=lambda r:r[0])
# group by material / t
groups=[]
for r in res:
    if groups and abs(groups[-1][-1][0]-r[0])<1e-3 and groups[-1][-1][2]==r[2]: groups[-1].append(r)
    else: groups.append([r])
tree=cKDTree(P)
out=[]
for g in groups:
    tt,c,mi,X=g[0]
    cand=np.unique(C[[r[1] for r in g]].ravel())
    dd=np.linalg.norm(P[cand]-X,axis=1); k=cand[np.argmin(dd)]
    # dist along borehole md from start (top?) 
    print(f'MatID {mi}: t={tt:.3f} m  X={X}  cells={[r[1] for r in g]}  nearest node {k} dist {dd.min():.3f}  node={P[k]}')
    out.append((mi,tt,X,k,dd.min(),[r[1] for r in g]))
# also distance of borehole line to each fracture overall
for mi in np.unique(mat):
    nodes=np.unique(C[mat==mi].ravel()); Q=P[nodes]
    dist=np.linalg.norm(np.cross(Q-A,u),axis=1); tq=(Q-A)@u
    print('Mat',mi,'nodes',len(nodes),'min dist to BH axis',dist.min().round(3),'at t',tq[np.argmin(dist)].round(2))
print('---')
for mi in np.unique(mat):
    nodes=np.unique(C[mat==mi].ravel()); Q=P[nodes]; c0=Q.mean(0)
    _,sv,vt=np.linalg.svd(Q-c0); n=vt[2]
    tp=((c0-A)@n)/(u@n); X=A+tp*u
    # distance of X to nearest fracture node
    print(f'Mat {mi}: planarity {sv[2]:.2e}, line-plane t={tp:.2f} m (BH length {L:.2f}), X to nearest frac node {np.linalg.norm(Q-X,axis=1).min():.2f} m')
for k in [270,271]:
    cells=np.where((C==k).any(1))[0]; print('node',k,'in MatIDs',np.unique(mat[cells]),'n cells',len(cells))
# which BH10 pts are these
print('BH10 pts idx for 270/271:', [int(np.argmin(np.linalg.norm(B-P[k],axis=1))) for k in [270,271]])
