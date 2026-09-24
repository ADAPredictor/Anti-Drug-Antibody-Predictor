"""Experimental-structure support: PDB fetching and solvent exposure.

Residue relative solvent-accessible area (RSA) is computed with a
grid-accelerated Shrake-Rupley algorithm (stdlib only).  Exposed residues
(RSA > 0.25) are used by the B-cell component because antibody-accessible
epitopes must lie on the protein surface.
"""
import math

from .config import (ALPHAFOLD_API_URL, RCSB_PDB_URL, STRUCTURE_DIR,
                     HTTP_TIMEOUT, USER_AGENT)

THREE_TO_ONE = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
    "MSE": "M",
}

# Empirical max SASA of residues in extended tripeptide (Tien et al. 2013)
MAX_ASA = {
    "A": 129.0, "R": 274.0, "N": 195.0, "D": 193.0, "C": 167.0,
    "Q": 225.0, "E": 223.0, "G": 104.0, "H": 224.0, "I": 182.0,
    "L": 180.0, "K": 236.0, "M": 200.0, "F": 212.0, "P": 136.0,
    "S": 155.0, "T": 172.0, "W": 258.0, "Y": 229.0, "V": 160.0,
}

ATOM_RADIUS = {"N": 1.65, "C": 1.76, "O": 1.40, "S": 1.85}
PROBE = 1.4
N_POINTS = 60
EXPOSED_CUTOFF = 0.25


def _sphere_points(n=N_POINTS):
    pts = []
    golden = math.pi * (3.0 - math.sqrt(5.0))
    for i in range(n):
        z = 1.0 - 2.0 * i / (n - 1)
        r = math.sqrt(max(0.0, 1.0 - z * z))
        theta = golden * i
        pts.append((math.cos(theta) * r, math.sin(theta) * r, z))
    return pts


def fetch_pdb_text(uniprot_acc=None, pdb_id=None):
    """Return (pdb_text, source_label). Network/cache errors raise IOError."""
    import requests

    headers = {"User-Agent": USER_AGENT, "Accept": "application/json"}
    if uniprot_acc:
        local = STRUCTURE_DIR / ("AF-%s.pdb" % uniprot_acc.upper())
        if local.exists():
            return local.read_text(encoding="utf-8"), "AlphaFold DB (cache)"
        resp = requests.get(ALPHAFOLD_API_URL.format(accession=uniprot_acc.upper()),
                            headers=headers, timeout=HTTP_TIMEOUT)
        resp.raise_for_status()
        entries = resp.json()
        pdb_url = entries[0]["pdbUrl"]
        pdb_resp = requests.get(pdb_url, headers={"User-Agent": USER_AGENT},
                                timeout=HTTP_TIMEOUT)
        pdb_resp.raise_for_status()
        local.write_text(pdb_resp.text, encoding="utf-8")
        return pdb_resp.text, "AlphaFold DB"
    if pdb_id:
        local = STRUCTURE_DIR / ("%s.pdb" % pdb_id.upper())
        if local.exists():
            return local.read_text(encoding="utf-8"), "RCSB PDB (cache)"
        resp = requests.get(RCSB_PDB_URL.format(pdb_id=pdb_id.upper()),
                            headers={"User-Agent": USER_AGENT}, timeout=HTTP_TIMEOUT)
        resp.raise_for_status()
        local.write_text(resp.text, encoding="utf-8")
        return resp.text, "RCSB PDB"
    raise ValueError("需要提供 uniprot_acc 或 pdb_id")


def parse_pdb(pdb_text):
    """Parse ATOM records -> ordered residues and atom coordinates.

    Returns: list of dicts(chain, resnum, icode, aa, atoms=[(el,x,y,z)...])
    """
    residues = []
    index = {}
    past_first_model = False
    for line in pdb_text.splitlines():
        if line.startswith("MODEL"):
            if past_first_model:
                break
            past_first_model = False  # first model starts
        if line.startswith("ENDMDL"):
            past_first_model = True   # ignore further (symmetry/copies)
            continue
        if not line.startswith("ATOM"):
            continue
        element = line[76:78].strip().upper() or line[12:16].strip()[0]
        if element not in ATOM_RADIUS:
            element = "C"
        chain = line[21].strip() or "A"
        resnum = int(line[22:26])
        icode = line[26].strip()
        resname = line[17:20].strip().upper()
        aa = THREE_TO_ONE.get(resname)
        if aa is None:
            continue
        key = (chain, resnum, icode)
        if key not in index:
            index[key] = len(residues)
            residues.append({"chain": chain, "resnum": resnum, "icode": icode,
                             "aa": aa, "atoms": []})
        x = float(line[30:38])
        y = float(line[38:46])
        z = float(line[46:54])
        residues[index[key]]["atoms"].append((element, x, y, z))
    return residues


def _shrake_rupley(residues):
    """Per-residue SASA in A^2."""
    atoms = []
    atom_res = []
    for ri, res in enumerate(residues):
        for el, x, y, z in res["atoms"]:
            atoms.append((x, y, z, ATOM_RADIUS[el] + PROBE))
            atom_res.append(ri)

    # spatial hashing
    cell_size = 2 * max(a[3] for a in atoms) if atoms else 1.0
    cells = {}
    for ai, (x, y, z, r) in enumerate(atoms):
        key = (int(math.floor(x / cell_size)),
               int(math.floor(y / cell_size)),
               int(math.floor(z / cell_size)))
        cells.setdefault(key, []).append(ai)

    pts = _sphere_points()
    sasa_atom = [0.0] * len(atoms)
    neighbor_offsets = [(dx, dy, dz) for dx in (-1, 0, 1)
                        for dy in (-1, 0, 1) for dz in (-1, 0, 1)]
    pi4 = 4.0 * math.pi

    for ai, (x, y, z, r) in enumerate(atoms):
        cx, cy, cz = (int(math.floor(x / cell_size)),
                      int(math.floor(y / cell_size)),
                      int(math.floor(z / cell_size)))
        # pre-unpack neighbours into flat lists for the hot loop
        nx_, ny_, nz_, nr2_ = [], [], [], []
        for dx, dy, dz in neighbor_offsets:
            for aj in cells.get((cx + dx, cy + dy, cz + dz), ()):
                if aj != ai:
                    ox, oy, oz, orad = atoms[aj]
                    nx_.append(ox)
                    ny_.append(oy)
                    nz_.append(oz)
                    nr2_.append(orad * orad)
        accessible = 0
        nnear = len(nx_)
        for px, py, pz in pts:
            sx, sy, sz = x + r * px, y + r * py, z + r * pz
            buried = False
            for k in range(nnear):
                ddx, ddy, ddz = sx - nx_[k], sy - ny_[k], sz - nz_[k]
                if ddx * ddx + ddy * ddy + ddz * ddz < nr2_[k]:
                    buried = True
                    break
            if not buried:
                accessible += 1
        sasa_atom[ai] = pi4 * r * r * accessible / N_POINTS

    sasa_res = [0.0] * len(residues)
    for ai, ri in enumerate(atom_res):
        sasa_res[ri] += sasa_atom[ai]
    return sasa_res


def compute_exposure(pdb_text, chains):
    """Align structure residues to query chains; return RSA/exposure track.

    Output rsa aligns to concatenated query sequence (None where unmapped).
    """
    residues = parse_pdb(pdb_text)
    if not residues:
        raise ValueError("PDB 中未解析到 ATOM 记录")

    # cheap per-structure-chain sequences (no SASA yet) for sequence matching
    seq_by_chain = {}
    for res in residues:
        seq_by_chain.setdefault(res["chain"], []).append(res["aa"])
    seq_by_chain = {k: "".join(v) for k, v in seq_by_chain.items()}

    def _match_chain(query_seq):
        for sc_id, sseq in seq_by_chain.items():
            pos = sseq.find(query_seq)
            if pos >= 0:
                return sc_id, pos
        for sc_id, sseq in seq_by_chain.items():
            if abs(len(sseq) - len(query_seq)) <= 3:
                return sc_id, 0
        return None

    # only calculate SASA on structure chains we actually need
    needed = set()
    for chain in chains:
        m = _match_chain(chain.seq)
        if m:
            needed.add(m[0])
    selected = [r for r in residues if r["chain"] in needed] or residues
    sasa = _shrake_rupley(selected)

    struct_chains = {}
    for ri, res in enumerate(selected):
        struct_chains.setdefault(res["chain"], []).append((res, sasa[ri]))

    total_len = sum(len(c.seq) for c in chains)
    rsa_track = [None] * total_len
    offset = 0
    mapped = 0
    exposure_flags = [0] * total_len
    for chain in chains:
        # find structure chain whose sequence contains this query chain
        best = None
        for sc_id, items in struct_chains.items():
            sseq = "".join(r["aa"] for r, _ in items)
            pos = sseq.find(chain.seq)
            if pos >= 0:
                best = (sc_id, pos)
                break
            # tolerate X by also trying longest chain if lengths match closely
            if best is None and abs(len(sseq) - len(chain.seq)) <= 3:
                best = (sc_id, 0)
        if best:
            sc_id, pos = best
            items = struct_chains[sc_id]
            for i, aa in enumerate(chain.seq):
                si = pos + i
                if 0 <= si < len(items):
                    res, rs = items[si]
                    if aa == "X" or res["aa"] == aa:
                        ref = MAX_ASA.get(aa, 200.0)
                        rsa = round(rs / ref, 3)
                        rsa_track[offset + i] = rsa
                        exposure_flags[offset + i] = 1 if rsa > EXPOSED_CUTOFF else 0
                        mapped += 1
        offset += len(chain.seq)

    mapped_fraction = round(mapped / total_len, 3) if total_len else 0.0
    mapped_rsas = [v for v in rsa_track if v is not None]
    return {
        "rsa_track": rsa_track,
        "exposed_track": exposure_flags,
        "mapped_fraction": mapped_fraction,
        "exposed_fraction": round(
            sum(v for v in exposure_flags) / len(mapped_rsas), 3
        ) if mapped_rsas else None,
        "n_residues_model": len(residues),
    }
