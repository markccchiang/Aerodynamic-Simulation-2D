const DEFAULTS = { alpha: 5, re_log: 6, m: 2, p: 4, t: 12 };
const PANELS = 160;
const flowCfg = { responsive: true, displayModeBar: false };
const cpCfg   = { responsive: true, displayModeBar: false };

// ---- read control values ----
function controls() {
    const o = {};
    document.querySelectorAll('.ctl').forEach(c => {
        o[c.dataset.key] = parseFloat(c.querySelector('input').value);
    });
    return o;
}
function fmt(k, v) {
    if (k === 'alpha') return v.toFixed(1) + '°';
    if (k === 're_log') return (10 ** v).toExponential(1).replace('e+', 'e');
    return String(v | 0);
}
function refreshLabels() {
    const c = controls();
    document.querySelectorAll('.ctl').forEach(el => {
        const k = el.dataset.key;
        el.querySelector('.val').textContent = fmt(k, c[k]);
    });
}

// ---- plotting ----
function drawFlow(d) {
    const sx = [], sy = [];
    for (const ln of d.streamlines) {
        for (let k = 0; k < ln.x.length; k++) { sx.push(ln.x[k]); sy.push(ln.y[k]); }
        sx.push(null); sy.push(null);
    }
    const heat = {
        type: 'heatmap', x: d.field.x, y: d.field.y, z: d.field.cp,
        colorscale: 'RdBu', zmin: -3, zmax: 1, zsmooth: 'best',
        colorbar: { title: { text: 'C<sub>p</sub>', side: 'right', font: { color: '#cdd6e6' } },
                    thickness: 12, len: 0.9, tickfont: { color: '#cdd6e6' } },
        hoverinfo: 'skip',
    };
    const streams = {
        type: 'scatter', mode: 'lines', x: sx, y: sy,
        line: { color: 'rgba(15,18,30,0.45)', width: 1 }, hoverinfo: 'skip',
    };
    const body = {
        type: 'scatter', mode: 'lines', x: d.geom.x, y: d.geom.y, fill: 'toself',
        fillcolor: '#e8eaed', line: { color: '#11141f', width: 1.4 }, hoverinfo: 'skip',
    };
    const layout = {
        margin: { l: 44, r: 10, t: 34, b: 38 },
        paper_bgcolor: '#1a2030', plot_bgcolor: '#1a2030',
        font: { color: '#cdd6e6' }, showlegend: false,
        title: { text: `${d.name}   α = ${d.alpha.toFixed(1)}°   Re = ${d.re.toExponential(1).replace('e+','e')}`,
                 font: { size: 14 } },
        xaxis: { title: 'x / c', range: [-0.6, 1.6], zeroline: false, gridcolor: '#283047' },
        yaxis: { title: 'y / c', range: [-0.7, 0.7], scaleanchor: 'x', scaleratio: 1,
                 zeroline: false, gridcolor: '#283047' },
    };
    Plotly.react('flow', [heat, streams, body], layout, flowCfg);
}

function drawCp(d) {
    const upper = { type: 'scatter', mode: 'lines', x: d.surface.x_upper, y: d.surface.cp_upper,
                    name: 'upper', line: { color: '#ff6b5e', width: 2 } };
    const lower = { type: 'scatter', mode: 'lines', x: d.surface.x_lower, y: d.surface.cp_lower,
                    name: 'lower', line: { color: '#4a9eff', width: 2 } };
    const shapes = [];
    const mark = (x, color) => { if (x != null) shapes.push({
        type: 'line', x0: x, x1: x, yref: 'paper', y0: 0, y1: 1,
        line: { color, width: 1.2, dash: 'dot' } }); };
    mark(d.coeffs.x_tr_upper, '#ff6b5e');
    mark(d.coeffs.x_tr_lower, '#4a9eff');
    const layout = {
        margin: { l: 46, r: 12, t: 34, b: 38 },
        paper_bgcolor: '#1a2030', plot_bgcolor: '#1a2030',
        font: { color: '#cdd6e6' }, shapes,
        title: { text: 'Surface pressure  (··· transition)', font: { size: 13 } },
        xaxis: { title: 'x / c', range: [-0.02, 1.02], gridcolor: '#283047', zeroline: false },
        yaxis: { title: 'C<sub>p</sub>', autorange: 'reversed', gridcolor: '#283047', zeroline: true,
                 zerolinecolor: '#445566' },
        legend: { x: 0.98, y: 0.04, xanchor: 'right', yanchor: 'bottom', font: { size: 11 } },
    };
    Plotly.react('cp', [upper, lower], layout, cpCfg);
}

function drawCoeffs(d) {
    const c = d.coeffs;
    document.getElementById('c-code').textContent = d.name;
    document.getElementById('c-cl').textContent  = c.cl.toFixed(3);
    document.getElementById('c-cd').textContent  = c.cd_visc.toFixed(4);
    document.getElementById('c-cm').textContent  = c.cm.toFixed(3);
    document.getElementById('c-cdp').textContent = c.cd_pressure.toFixed(5);
    const sep = document.getElementById('c-sep');
    if (c.separated) { sep.className = 'badge sep'; sep.textContent = 'trailing-edge separation'; }
    else { sep.className = 'badge attached'; sep.textContent = 'attached'; }
}

// ---- airfoil source: NACA sliders, a bundled sample, or an upload ----
let source = { type: 'naca' };
const nacaKeys = ['m', 'p', 't'];

function setNacaEnabled(on) {
    document.querySelectorAll('.ctl').forEach(el => {
        if (nacaKeys.includes(el.dataset.key))
            el.classList.toggle('disabled', !on);
    });
}

async function loadSamples() {
    try {
        const res = await fetch('/api/samples');
        const { samples } = await res.json();
        const sel = document.getElementById('source');
        for (const s of samples) {
            const o = document.createElement('option');
            o.value = 'sample:' + s.key; o.textContent = s.name;
            sel.appendChild(o);
        }
    } catch (e) { console.error(e); }
}

document.getElementById('source').addEventListener('change', (e) => {
    const v = e.target.value;
    if (v === 'naca') { source = { type: 'naca' }; setNacaEnabled(true); }
    else if (v.startsWith('sample:')) {
        source = { type: 'sample', sample: v.slice(7), name: e.target.selectedOptions[0].textContent };
        setNacaEnabled(false);
    }
    document.getElementById('loaded').textContent = '';
    schedule();
});

document.getElementById('upload').addEventListener('change', (e) => {
    const file = e.target.files[0];
    if (!file) return;
    const reader = new FileReader();
    reader.onload = () => {
        source = { type: 'custom', dat: reader.result, name: file.name.replace(/\.[^.]+$/, '') };
        // Add/replace an "uploaded" option and select it.
        const sel = document.getElementById('source');
        let opt = sel.querySelector('option[value="upload"]');
        if (!opt) { opt = document.createElement('option'); opt.value = 'upload'; sel.appendChild(opt); }
        opt.textContent = '↑ ' + file.name;
        sel.value = 'upload';
        setNacaEnabled(false);
        schedule();
    };
    reader.readAsText(file);
    e.target.value = '';  // allow re-uploading the same file
});

// ---- fetch + render (debounced, drops stale responses) ----
let reqId = 0, timer = null;
async function solve() {
    const c = controls();
    const mine = ++reqId;
    try {
        let res;
        if (source.type === 'naca') {
            res = await fetch('/api/solve?' + new URLSearchParams({ ...c, panels: PANELS }));
        } else {
            const body = { alpha: c.alpha, re_log: c.re_log, panels: PANELS, name: source.name };
            if (source.type === 'sample') body.sample = source.sample; else body.dat = source.dat;
            res = await fetch('/api/solve_custom', {
                method: 'POST', headers: { 'Content-Type': 'application/json' }, body: JSON.stringify(body),
            });
        }
        if (!res.ok) {
            const msg = await res.json().catch(() => ({}));
            throw new Error(msg.detail || ('HTTP ' + res.status));
        }
        const d = await res.json();
        if (mine !== reqId) return;            // a newer request superseded this one
        if (source.type === 'naca') {
            // Mirror the server's camber-position rule (P snaps to 1 when M>0, P=0).
            const pIn = document.querySelector('.ctl[data-key="p"] input');
            if (parseInt(pIn.value, 10) !== d.p) { pIn.value = d.p; refreshLabels(); }
        } else {
            document.getElementById('loaded').textContent =
                `loaded ${d.n_points} points → ${d.geom.x.length - 1} panels`;
        }
        drawFlow(d); drawCp(d); drawCoeffs(d);
    } catch (e) {
        document.getElementById('loaded').textContent = '⚠ ' + e.message;
        console.error(e);
    }
}
function schedule() { refreshLabels(); clearTimeout(timer); timer = setTimeout(solve, 60); }

document.querySelectorAll('.ctl input').forEach(inp => inp.addEventListener('input', schedule));
document.getElementById('reset').addEventListener('click', () => {
    document.querySelectorAll('.ctl').forEach(el => {
        el.querySelector('input').value = DEFAULTS[el.dataset.key];
    });
    source = { type: 'naca' };
    document.getElementById('source').value = 'naca';
    document.getElementById('loaded').textContent = '';
    setNacaEnabled(true);
    schedule();
});

loadSamples();
refreshLabels();
solve();
