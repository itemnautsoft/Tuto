import struct, re
class Blend:
    def __init__(s, path):
        d = s.d = open(path,'rb').read(); assert d[:7]==b'BLENDER'
        s.blocks=[]; o=12; s.byptr={}
        while o < len(d):
            code=d[o:o+4]; size,=struct.unpack_from('<i',d,o+4); ptr,=struct.unpack_from('<Q',d,o+8); sdna,cnt=struct.unpack_from('<ii',d,o+16)
            b=dict(code=code.rstrip(b'\0').decode('latin1'),size=size,ptr=ptr,sdna=sdna,cnt=cnt,off=o+24); s.blocks.append(b)
            if ptr: s.byptr[ptr]=b
            if code==b'ENDB': break
            o+=24+size
        dna=[b for b in s.blocks if b['code']=='DNA1'][0]; p=dna['off']+8
        def al(p): return (p+3)&~3
        n,=struct.unpack_from('<i',d,p); p+=4; s.names=[]
        for i in range(n): e=d.index(b'\0',p); s.names.append(d[p:e].decode()); p=e+1
        p=al(p)+4; n,=struct.unpack_from('<i',d,p); p+=4; s.types=[]
        for i in range(n): e=d.index(b'\0',p); s.types.append(d[p:e].decode()); p=e+1
        p=al(p)+4; s.tlen=struct.unpack_from('<%dh'%n,d,p); p+=2*n
        p=al(p)+4; n,=struct.unpack_from('<i',d,p); p+=4; s.structs=[]; s.sidx={}
        for i in range(n):
            t,nf=struct.unpack_from('<hh',d,p); p+=4; f=[]; off=0
            for j in range(nf):
                ft,fn=struct.unpack_from('<hh',d,p); p+=4; name=s.names[fn]
                ptr=name.startswith('*') or name.startswith('(*'); arr=1
                for a in re.findall(r'\[(\d+)\]',name): arr*=int(a)
                sz=(8 if ptr else s.tlen[ft])*arr; base=re.sub(r'[\*\(\)]|\[.*','',name)
                f.append((base,s.types[ft],ptr,arr,off,sz,name)); off+=sz
            s.structs.append((s.types[t],f)); s.sidx[s.types[t]]=i
    def fields(s, tname): return {f[0]:f for f in s.structs[s.sidx[tname]][1]}
    def get(s, off, tname, field):
        f=s.fields(tname)[field]; base,ft,ptr,arr,fo,sz,_=f; o=off+fo; d=s.d
        if ptr: return struct.unpack_from('<Q',d,o)[0] if arr==1 else struct.unpack_from('<%dQ'%arr,d,o)
        fmt={'int':'i','float':'f','short':'h','char':'b','uchar':'B','int64_t':'q','uint64_t':'Q','double':'d','int8_t':'b','uint8_t':'B','ushort':'H','uint':'I','int16_t':'h','int32_t':'i','uint32_t':'I'}.get(ft)
        if ft=='char' and arr>1: return d[o:o+arr].split(b'\0')[0].decode('latin1')
        if fmt: v=struct.unpack_from('<%d%s'%(arr,fmt),d,o); return v[0] if arr==1 else v
        return o
    def blk(s, ptr): return s.byptr.get(ptr)
