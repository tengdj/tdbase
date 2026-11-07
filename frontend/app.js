const $ = (id) => document.getElementById(id);
const form = $('query-form');
const query = $('query');
const dimensionInputs = [...document.querySelectorAll('input[name="dimension"]')];
const backendInputs = [...document.querySelectorAll('input[name="backend"]')];
let lastResult = null;

const queryOptions = {
  '2d': [['intersect', '相交 Intersect'], ['intersection', '交集面积 Intersection'], ['within', '距离 DWithin']],
  '3d': [['intersect', '相交 Intersect'], ['within', '距离 DWithin'], ['nn', '近邻 KNN']],
};

function selected(name) { return document.querySelector(`input[name="${name}"]:checked`).value; }
function setVisible(id, visible) { $(id).classList.toggle('hidden', !visible); }

function refreshForm() {
  const dimension = selected('dimension');
  const current = query.value;
  query.replaceChildren(...queryOptions[dimension].map(([value, label]) => new Option(label, value)));
  query.value = queryOptions[dimension].some(([value]) => value === current) ? current : queryOptions[dimension][0][0];
  const is3d = dimension === '3d';
  const isWithin = query.value === 'within';
  const isKnn = is3d && query.value === 'nn';
  $('run-label').textContent = `运行${is3d ? '三维' : '二维'}查询`;
  $('distance-unit').textContent = is3d ? '三维单位' : 'km';
  setVisible('k-row', isKnn);
  setVisible('distance-field', isWithin);
  setVisible('vpr-group', !is3d);
  setVisible('bitwidth-group', !is3d);
  setVisible('compute-threads-group', is3d);
  setVisible('target-object-group', is3d);
  $('target').placeholder = is3d ? '…/data/target3d.dt' : '…/data/target2d.dt';
  document.querySelectorAll('.segmented').forEach(group => group.querySelectorAll('.seg').forEach(label => label.classList.toggle('active', label.querySelector('input').checked)));
}

dimensionInputs.forEach(input => input.addEventListener('change', refreshForm));
backendInputs.forEach(input => input.addEventListener('change', refreshForm));
query.addEventListener('change', refreshForm);
refreshForm();

function showError(message) {
  $('empty-state').classList.remove('hidden');
  $('result-state').classList.add('hidden');
  $('empty-state').innerHTML = `<div class="empty-icon error">!</div><h3>查询未完成</h3><p>${escapeHtml(message)}</p>`;
}
function escapeHtml(value) { return String(value).replace(/[&<>'"]/g, ch => ({'&':'&amp;','<':'&lt;','>':'&gt;',"'":'&#39;','"':'&quot;'}[ch])); }
function showResult(result) {
  lastResult = result;
  $('empty-state').classList.add('hidden');
  $('result-state').classList.remove('hidden');
  $('result-status').textContent = result.ok ? '完成' : `失败 (${result.exit_code})`;
  $('result-status').className = result.ok ? 'success' : 'failure';
  $('result-count').textContent = result.count == null ? '—' : result.count;
  $('result-duration').textContent = `${result.duration_seconds}s`;
  $('result-command').textContent = result.command.map(part => /\s/.test(part) ? JSON.stringify(part) : part).join(' ');
  $('result-output').textContent = result.stdout || result.stderr || '没有输出';
  document.querySelectorAll('.tab').forEach(tab => tab.classList.toggle('active', tab.dataset.output === 'stdout'));
}
document.querySelectorAll('.tab').forEach(tab => tab.addEventListener('click', () => {
  document.querySelectorAll('.tab').forEach(item => item.classList.toggle('active', item === tab));
  $('result-output').textContent = lastResult ? (lastResult[tab.dataset.output] || '没有输出') : '';
}));

form.addEventListener('submit', async event => {
  event.preventDefault();
  const button = form.querySelector('.run-button');
  button.disabled = true;
  button.classList.add('running');
  $('run-label').textContent = '正在运行…';
  try {
    const payload = {
      dimension: selected('dimension'), backend: selected('backend'), query: query.value,
      source: $('source').value, target: $('target').value, distance: $('distance').value,
      target_object: $('target-object').value,
      k: $('k').value, lods: $('lods').value, threads: $('threads').value,
      vpr: $('vpr').value, bitwidth: $('bitwidth').value, 'compute_threads': $('compute-threads').value,
      hierarchy: $('hierarchy').checked, approximation: $('approximation').checked,
      print_result: $('print-result').checked,
    };
    const response = await fetch('/api/query', {method: 'POST', headers: {'Content-Type': 'application/json'}, body: JSON.stringify(payload)});
    const result = await response.json();
    if (!response.ok && result.error) throw new Error(result.error);
    showResult(result);
  } catch (error) { showError(error.message); }
  finally { button.disabled = false; button.classList.remove('running'); refreshForm(); }
});

async function health() {
  try {
    const response = await fetch('/api/health');
    const data = await response.json();
    $('health').innerHTML = `<span class="health-dot ${data.gpu_available ? 'online' : 'warn'}"></span><span>CPU ${data.cpu_available ? '就绪' : '缺失'} · GPU ${data.gpu_available ? '就绪' : '未编译'}</span>`;
  } catch (error) { $('health').innerHTML = '<span class="health-dot warn"></span><span>服务连接失败</span>'; }
}
health();
