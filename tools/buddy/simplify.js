// Fast quadric mesh simplification (after Sven Forstmann), with a "class" lock: edges whose two ends have different colour classes never collapse
const fs = require('fs');
const [,, inF, outF, targetS] = process.argv; const target = +targetS;
const buf = fs.readFileSync(inF); const nv = buf.readUInt32LE(0), nt = buf.readUInt32LE(4);
const P0 = new Float32Array(buf.buffer.slice(buf.byteOffset + 8, buf.byteOffset + 8 + nv * 12));
const T0 = new Int32Array(buf.buffer.slice(buf.byteOffset + 8 + nv * 12, buf.byteOffset + 8 + nv * 12 + nt * 12));
const C0 = new Uint8Array(buf.buffer.slice(buf.byteOffset + 8 + nv * 12 + nt * 12, buf.byteOffset + 8 + nv * 12 + nt * 12 + nv));
let V = []; for (let i = 0; i < nv; i++) V.push({ p: [P0[i * 3], P0[i * 3 + 1], P0[i * 3 + 2]], q: new Float64Array(10), ts: 0, tc: 0, border: 0, c: C0[i] });
let T = []; for (let i = 0; i < nt; i++) T.push({ v: [T0[i * 3], T0[i * 3 + 1], T0[i * 3 + 2]], err: [0, 0, 0, 0], del: 0, dirty: 0, n: [0, 0, 0] });
let R = [];   // refs {t, tv}
const sub = (a, b) => [a[0] - b[0], a[1] - b[1], a[2] - b[2]], cross = (a, b) => [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]], dot = (a, b) => a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
const norm = a => { const l = Math.hypot(a[0], a[1], a[2]) || 1; return [a[0] / l, a[1] / l, a[2] / l]; };
const planeQ = (a, b, c, d) => new Float64Array([a * a, a * b, a * c, a * d, b * b, b * c, b * d, c * c, c * d, d * d]);
const det = (q, a11, a12, a13, a21, a22, a23, a31, a32, a33) => q[a11] * q[a22] * q[a33] + q[a13] * q[a21] * q[a32] + q[a12] * q[a23] * q[a31] - q[a13] * q[a22] * q[a31] - q[a11] * q[a23] * q[a32] - q[a12] * q[a21] * q[a33];
const vErr = (q, x, y, z) => q[0] * x * x + 2 * q[1] * x * y + 2 * q[2] * x * z + 2 * q[3] * x + q[4] * y * y + 2 * q[5] * y * z + 2 * q[6] * y + q[7] * z * z + 2 * q[8] * z + q[9];
function calcErr(i0, i1, out) {
  const q = new Float64Array(10); for (let k = 0; k < 10; k++) q[k] = V[i0].q[k] + V[i1].q[k];
  const border = V[i0].border & V[i1].border; let e;
  const d = det(q, 0, 1, 2, 1, 4, 5, 2, 5, 7);
  if (d !== 0 && !border) {
    out[0] = -1 / d * det(q, 1, 2, 3, 4, 5, 6, 5, 7, 8); out[1] = 1 / d * det(q, 0, 2, 3, 1, 5, 6, 2, 7, 8); out[2] = -1 / d * det(q, 0, 1, 3, 1, 4, 6, 2, 5, 8);
    e = vErr(q, out[0], out[1], out[2]);
  } else {
    const p1 = V[i0].p, p2 = V[i1].p, p3 = [(p1[0] + p2[0]) / 2, (p1[1] + p2[1]) / 2, (p1[2] + p2[2]) / 2];
    const e1 = vErr(q, ...p1), e2 = vErr(q, ...p2), e3 = vErr(q, ...p3); e = Math.min(e1, e2, e3);
    const p = e1 === e ? p1 : e2 === e ? p2 : p3; out[0] = p[0]; out[1] = p[1]; out[2] = p[2];
  }
  return e;
}
function flipped(p, i0, i1, v0, deleted) {
  for (let k = 0; k < v0.tc; k++) {
    const r = R[v0.ts + k], t = T[r.t]; if (t.del) continue;
    const s = r.tv, id1 = t.v[(s + 1) % 3], id2 = t.v[(s + 2) % 3];
    if (id1 === i1 || id2 === i1) { deleted[k] = 1; continue; }
    const d1 = norm(sub(V[id1].p, p)), d2 = norm(sub(V[id2].p, p));
    if (Math.abs(dot(d1, d2)) > 0.999) return true;
    const n = norm(cross(d1, d2)); deleted[k] = 0;
    if (dot(n, t.n) < 0.2) return true;
  }
  return false;
}
let delTris = 0;
function updateTris(i0, v, deleted) {
  const p = [0, 0, 0];
  for (let k = 0; k < v.tc; k++) {
    const r = R[v.ts + k], t = T[r.t]; if (t.del) continue;
    if (deleted[k]) { t.del = 1; delTris++; continue; }
    t.v[r.tv] = i0; t.dirty = 1;
    t.err[0] = calcErr(t.v[0], t.v[1], p); t.err[1] = calcErr(t.v[1], t.v[2], p); t.err[2] = calcErr(t.v[2], t.v[0], p); t.err[3] = Math.min(t.err[0], t.err[1], t.err[2]);
    R.push({ t: r.t, tv: r.tv });
  }
}
function updateMesh(iter) {
  if (iter > 0) { T = T.filter(t => !t.del); }
  if (iter === 0) {
    for (const v of V) v.q.fill(0);
    for (const t of T) {
      const p = t.v.map(i => V[i].p), n = norm(cross(sub(p[1], p[0]), sub(p[2], p[0]))); t.n = n;
      const q = planeQ(n[0], n[1], n[2], -dot(n, p[0])); for (const i of t.v) for (let k = 0; k < 10; k++) V[i].q[k] += q[k];
    }
    const pp = [0, 0, 0];
    for (const t of T) { for (let j = 0; j < 3; j++) t.err[j] = calcErr(t.v[j], t.v[(j + 1) % 3], pp); t.err[3] = Math.min(t.err[0], t.err[1], t.err[2]); }
  }
  for (const v of V) { v.ts = 0; v.tc = 0; }
  for (const t of T) for (const i of t.v) V[i].tc++;
  let ts = 0; for (const v of V) { v.ts = ts; ts += v.tc; v.tc = 0; }
  R = new Array(T.length * 3);
  T.forEach((t, ti) => { for (let j = 0; j < 3; j++) { const v = V[t.v[j]]; R[v.ts + v.tc] = { t: ti, tv: j }; v.tc++; } });
  if (iter === 0) {
    for (let vi = 0; vi < V.length; vi++) {
      const v = V[vi], cnt = new Map();
      for (let k = 0; k < v.tc; k++) { const t = T[R[v.ts + k].t]; for (const id of t.v) cnt.set(id, (cnt.get(id) || 0) + 1); }
      for (const [id, c] of cnt) if (c === 1) V[id].border = 1;
    }
  }
}
let triCount = T.length;
for (let iter = 0; iter < 200; iter++) {
  if (triCount - delTris <= target) break;
  if (iter % 5 === 0) { updateMesh(iter); }
  for (const t of T) t.dirty = 0;
  const thr = 1e-9 * Math.pow(iter + 3, 7);
  for (const t of T) {
    if (t.err[3] > thr || t.del || t.dirty) continue;
    for (let j = 0; j < 3; j++) {
      if (t.err[j] >= thr) continue;
      const i0 = t.v[j], i1 = t.v[(j + 1) % 3], v0 = V[i0], v1 = V[i1];
      if (v0.border || v1.border) continue;
      if (v0.c !== v1.c) continue;
      const p = [0, 0, 0]; calcErr(i0, i1, p);
      const d0 = new Uint8Array(v0.tc), d1 = new Uint8Array(v1.tc);
      if (flipped(p, i0, i1, v0, d0) || flipped(p, i1, i0, v1, d1)) continue;
      v0.p = p; for (let k = 0; k < 10; k++) v0.q[k] += v1.q[k];
      const ts = R.length; updateTris(i0, v0, d0); updateTris(i0, v1, d1);
      const tc = R.length - ts;
      if (tc <= v0.tc) { for (let k = 0; k < tc; k++) R[v0.ts + k] = R[ts + k]; R.length = ts; }
      else v0.ts = ts;
      v0.tc = tc; break;
    }
    if (triCount - delTris <= target) break;
  }
  if (iter % 5 === 4) process.stderr.write(`iter ${iter} tris ${triCount - delTris}\n`);
}
// compact
T = T.filter(t => !t.del); const map = new Int32Array(V.length).fill(-1); const outV = [];
for (const t of T) for (let j = 0; j < 3; j++) { const i = t.v[j]; if (map[i] < 0) { map[i] = outV.length; outV.push(i); } t.v[j] = map[i]; }
const ob = Buffer.alloc(8 + outV.length * 13 + T.length * 12); ob.writeUInt32LE(outV.length, 0); ob.writeUInt32LE(T.length, 4); let o = 8;
for (const i of outV) { for (let k = 0; k < 3; k++) { ob.writeFloatLE(V[i].p[k], o); o += 4; } }
for (const t of T) for (let j = 0; j < 3; j++) { ob.writeInt32LE(t.v[j], o); o += 4; }
for (const i of outV) { ob.writeUInt8(V[i].c, o); o++; }
fs.writeFileSync(outF, ob); console.log('out verts', outV.length, 'tris', T.length);
