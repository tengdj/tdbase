const $ = (id) => document.getElementById(id);
const form = $('query-form');
const query = $('query');
const dimensionInputs = [...document.querySelectorAll('input[name="dimension"]')];
const backendInputs = [...document.querySelectorAll('input[name="backend"]')];
let lastResult = null;
let visualizationCleanup = null;
let visualizationKey = null;
const distanceValues = { '2d': '10', '3d': '100' };
let distanceDimension = null;

const queryOptions = {
  '2d': [['intersect', '相交 Intersect'], ['intersection', '交集面积 Intersection'], ['within', '距离 DWithin']],
  '3d': [['intersect', '相交 Intersect'], ['within', '距离 DWithin'], ['nn', '近邻 KNN']],
};

function projectPath(path) {
  const root = document.documentElement.dataset.projectRoot;
  return root ? `${root.replace(/\/$/, '')}/${path}` : path;
}

const threeDDataFiles = [
  'data/foo_n_nv50_nu200_vs100_r30.dt',
  'data/foo_v_nv50_nu200_vs100_r30.dt',
  'data/decoded_n_nv50_nu200_vs100_r30.dt',
  'data/decoded_v_nv50_nu200_vs100_r30.dt',
].map(projectPath);
const twoDDataFiles = [
  'data/complex.dt',
  'data/valid_lakes_polygons_100.dt',
  'data/oligo1_prj_div.dt',
  'data/oligo2_prj_div.dt',
  'data/lakes.dt',
].map(projectPath);
const twoDQueryDefaults = {
  intersect: { source: twoDDataFiles[0], target: twoDDataFiles[1] },
  intersection: { source: twoDDataFiles[2], target: twoDDataFiles[3] },
  within: { source: twoDDataFiles[4], target: twoDDataFiles[4] },
};
const dataPaths = {
  ...Object.fromEntries(Object.entries(twoDQueryDefaults).map(([type, paths]) => [`2d:${type}`, { ...paths }])),
  '3d': { source: threeDDataFiles[2], target: threeDDataFiles[2] },
};
const customPaths = {};
let pathContext = null;
const datasetStatsRequests = new Map();
const datasetStatsKeys = { source: null, target: null };
const datasetStatsTimers = {};

function datasetStatsKey(name) {
  const path = $(name).value.trim() || (name === 'target' ? $('source').value.trim() : '');
  return JSON.stringify([selected('dimension'), path]);
}

async function refreshDatasetStats(name) {
  const key = datasetStatsKey(name);
  if (datasetStatsKeys[name] === key) return;
  datasetStatsKeys[name] = key;
  const [dimension, path] = JSON.parse(key);
  const line = $(name + '-stats');
  line.dataset.state = path ? 'loading' : 'empty';
  line.textContent = path ? '正在读取数据统计…' : '选择文件后显示数据规模';
  line.removeAttribute('title');
  if (!path) return;
  if (!datasetStatsRequests.has(key)) {
    const request = fetch('/api/dataset-stats', {
      method: 'POST', headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify({ dimension, path }),
    }).then(async response => {
      const result = await response.json();
      if (!response.ok || !result.ok) throw new Error(result.error || '读取失败');
      return result;
    });
    datasetStatsRequests.set(key, request);
    // Keep in-flight requests shared; fresh selections revalidate file changes
    // against the server's size/mtime cache rather than stale browser counts.
    request.finally(() => datasetStatsRequests.delete(key)).catch(() => {});
  }
  try {
    const result = await datasetStatsRequests.get(key);
    if (datasetStatsKeys[name] !== key || datasetStatsKey(name) !== key) return;
    const number = value => new Intl.NumberFormat('zh-CN').format(value);
    // The totals are useful for validation, but they are hard to compare
    // between datasets with different object counts. Show the per-object
    // average instead. FARM's polygon vertex and edge counts are identical,
    // so the 2D row intentionally keeps only the average edge count.
    const average = (total) => result.objects > 0
      ? new Intl.NumberFormat('zh-CN', { minimumFractionDigits: 1, maximumFractionDigits: 1 }).format(total / result.objects)
      : '—';
    const values = dimension === '3d'
      ? [['物体', result.objects], ['平均顶点/物体', average(result.vertices)], ['平均面/物体', average(result.faces)], ['面总数', result.faces]]
      : [['多边形', result.objects], ['平均边/物体', average(result.edges)], ['顶点总数', result.vertices]];
    line.replaceChildren(...values.map(([label, value]) => {
      const item = document.createElement('span');
      const count = document.createElement('strong');
      count.textContent = typeof value === 'number' ? number(value) : value;
      item.append(document.createTextNode(label + ' '), count);
      return item;
    }));
    line.dataset.state = 'ready';
    line.title = dimension === '3d' ? '物体数和三角面保留总量；顶点和三角面同时显示每个物体的平均数量，统计完整精度（LOD 100）。'
      : '多边形数和顶点保留总量；二维顶点数与边数相同，因此平均值仅显示每个物体的平均边数。';
  } catch (error) {
    if (datasetStatsKeys[name] !== key || datasetStatsKey(name) !== key) return;
    line.dataset.state = 'error';
    line.textContent = '无法读取数据统计：' + error.message;
    line.title = error.message;
    // Allow a later selection/blur to retry transient failures.
    datasetStatsKeys[name] = null;
  }
}

function refreshAllDatasetStats() {
  for (const name of ['source', 'target']) refreshDatasetStats(name);
}

function filesForDimension(dimension) {
  return dimension === '3d' ? threeDDataFiles : twoDDataFiles;
}

function refreshPathChoice(name) {
  const custom = $(name + '-preset').value === 'custom';
  setVisible(name + '-preset-group', true);
  setVisible(name + '-path', custom);
  $(name + '-label').htmlFor = custom ? name : name + '-preset';
}

function refreshPresetOptions(name, dimension, custom) {
  const preset = $(name + '-preset');
  const current = $(name).value;
  const files = filesForDimension(dimension);
  preset.replaceChildren(
    ...files.map(path => new Option(path, path)),
    new Option('自定义路径…', 'custom'),
  );
  preset.value = !custom && files.includes(current) ? current : 'custom';
}

for (const name of ['source', 'target']) {
  const preset = $(name + '-preset');
  preset.addEventListener('change', () => {
    if (preset.value !== 'custom') $(name).value = preset.value;
    refreshPathChoice(name);
    if (preset.value === 'custom') $(name).focus();
    refreshAllDatasetStats();
  });
  $(name).addEventListener('input', () => {
    clearTimeout(datasetStatsTimers[name]);
    // Invalidate the old display immediately while a custom path is edited.
    datasetStatsKeys[name] = null;
    $(name + '-stats').dataset.state = 'loading';
    $(name + '-stats').textContent = '等待路径输入完成…';
    datasetStatsTimers[name] = setTimeout(refreshAllDatasetStats, 400);
  });
  $(name).addEventListener('change', refreshAllDatasetStats);
  $(name).addEventListener('blur', refreshAllDatasetStats);
}

function refreshDataPaths(dimension, queryValue) {
  const context = dimension === '2d' ? `2d:${queryValue}` : '3d';
  if (context !== pathContext) {
    if (pathContext) {
      customPaths[pathContext] = {};
      for (const name of ['source', 'target']) {
        dataPaths[pathContext][name] = $(name).value;
        customPaths[pathContext][name] = $(name + '-preset').value === 'custom';
      }
    }
    for (const name of ['source', 'target']) {
      $(name).value = dataPaths[context][name];
      $(name).placeholder = projectPath(`data/${name}${dimension}.dt`);
      refreshPresetOptions(name, dimension, customPaths[context]?.[name]);
      refreshPathChoice(name);
    }
    pathContext = context;
  }
}

function selected(name) { return document.querySelector(`input[name="${name}"]:checked`).value; }
function setVisible(id, visible) { $(id).classList.toggle('hidden', !visible); }

function refreshForm() {
  const dimension = selected('dimension');
  if (dimension !== distanceDimension) {
    if (distanceDimension) distanceValues[distanceDimension] = $('distance').value;
    $('distance').value = distanceValues[dimension];
    distanceDimension = dimension;
  }
  const current = query.value;
  query.replaceChildren(...queryOptions[dimension].map(([value, label]) => new Option(label, value)));
  query.value = queryOptions[dimension].some(([value]) => value === current) ? current : queryOptions[dimension][0][0];
  const is3d = dimension === '3d';
  const isWithin = query.value === 'within';
  const isKnn = is3d && query.value === 'nn';
  $('run-label').textContent = '运行查询';
  $('distance-unit').textContent = is3d ? '三维单位' : 'km';
  setVisible('k-field', isKnn);
  setVisible('distance-field', isWithin);
  setVisible('lod-row', is3d);
  setVisible('vpr-group', !is3d);
  setVisible('bitwidth-group', !is3d);
  setVisible('farm-options', !is3d);
  setVisible('compute-threads-group', is3d);
  setVisible('target-object-group', true);
  refreshDataPaths(dimension, query.value);
  refreshAllDatasetStats();
  document.querySelectorAll('.segmented').forEach(group => group.querySelectorAll('.seg').forEach(label => label.classList.toggle('active', label.querySelector('input').checked)));
}

dimensionInputs.forEach(input => input.addEventListener('change', refreshForm));
backendInputs.forEach(input => input.addEventListener('change', refreshForm));
query.addEventListener('change', refreshForm);
$('progressive-filter').addEventListener('change', refreshForm);
refreshForm();

function showError(message) {
  visualizationCleanup?.();
  visualizationCleanup = null;
  visualizationKey = null;
  lastResult = null;
  $('visualize-button').disabled = true;
  $('empty-state').classList.remove('hidden');
  $('result-state').classList.add('hidden');
  $('empty-state').innerHTML = `<div class="empty-icon error">!</div><h3>查询未完成</h3><p>${escapeHtml(message)}</p>`;
}
function escapeHtml(value) { return String(value).replace(/[&<>'"]/g, ch => ({'&':'&amp;','<':'&lt;','>':'&gt;',"'":'&#39;','"':'&quot;'}[ch])); }
function isNumberToken(value) { return value !== '' && Number.isFinite(Number(value)); }
function formatThreeDStats(line) {
  const values = line.trim().split(/\s+/);
  if (values.length !== 17 || !values.every(isNumberToken)) return null;
  const lines = ['【三维查询统计】', '每组数据依次为：候选对象数量、体素对数量'];
  for (let i = 0; i < 6; i += 1) {
    const lod = i === 0 ? '初始层' : `LOD ${i * 20}`;
    lines.push(`${lod}：候选对象 ${values[i * 2]}，体素对 ${values[i * 2 + 1]}`);
  }
  const labels = ['索引', '解码', '数据准备', '几何计算', '候选更新'];
  lines.push('阶段耗时（毫秒）：');
  lines.push(labels.map((label, i) => `${label} ${values[12 + i]} ms`).join('  ·  '));
  return lines.join('\n');
}
function formatOutput(result, outputName) {
  const raw = result[outputName] || '';
  if (!raw) return '没有输出';
  if (outputName === 'stdout' && result.dimension === '3d') {
    return raw.split(/\r?\n/).map(line => formatThreeDStats(line) || line).join('\n');
  }
  return raw;
}
function outputDescription(result, outputName) {
  if (outputName === 'stderr') {
    return '运行日志：记录数据加载、空间索引、LOD 解码、数据准备、几何计算、候选更新和最终统计。日志中的时间通常以毫秒显示，超过 1000 毫秒时会换算为秒。';
  }
  if (result.dimension === '3d') {
    const preview = result.output_omitted_pairs ? `当前显示前100对和后100对，中间省略 ${result.output_omitted_pairs.toLocaleString()} 对。` : '当前结果对数不超过200，已完整显示。';
    return `三维结果：普通行是“源对象 ID → 目标对象 ID”；末尾统计按 LOD 展示候选对象数、体素对数和各阶段耗时。${preview}`;
  }
  if (result.query === 'intersection') {
    const preview = result.output_omitted_pairs ? `当前显示前100对和后100对，中间省略 ${result.output_omitted_pairs.toLocaleString()} 对。` : '当前结果对数不超过200，已完整显示。';
    return `二维交集面积结果：每行依次为“源对象 ID、目标对象 ID、交集面积”；面积单位是输入坐标单位的平方。${preview}`;
  }
  const preview = result.output_omitted_pairs ? `当前显示前100对和后100对，中间省略 ${result.output_omitted_pairs.toLocaleString()} 对。` : '当前结果对数不超过200，已完整显示。';
  return `二维结果：每行依次为“源对象 ID、目标对象 ID”；运行日志和匹配数量请切换到“运行日志”查看。${preview}`;
}
function updateOutput(result, outputName) {
  $('output-help').textContent = outputDescription(result, outputName);
  $('result-output').textContent = formatOutput(result, outputName);
}
function formatQueryDuration(result, stage) {
  const seconds = result[`${stage}_seconds`];
  if (typeof seconds === 'number' && Number.isFinite(seconds) && seconds >= 0) {
    return seconds >= 1 ? `${seconds.toFixed(3)}s` : `${(seconds * 1000).toFixed(3)}ms`;
  }
  const logs = `${result.stderr || ''}\n${result.stdout || ''}`;
  const pattern = new RegExp(`^\\s*${stage}:\\s*([0-9]+(?:\\.[0-9]+)?(?:e[+-]?[0-9]+)?)\\s*(ms|s)\\b`, 'gim');
  const matches = [...logs.matchAll(pattern)];
  if (!matches.length) return '—';
  const match = matches[matches.length - 1];
  return `${match[1]}${match[2].toLowerCase()}`;
}
function updateTiming(result) {
  const summary = document.querySelector('.timing-summary');
  if (summary) {
    summary.classList.toggle('timing-3d', result.dimension === '3d');
    summary.classList.toggle('timing-2d', result.dimension === '2d');
  }
  const stages = ['preprocess', 'index', 'compute', 'rasterization'];
  const values = stages.map(stage => result[`${stage}_seconds`]).filter(value =>
    typeof value === 'number' && Number.isFinite(value) && value >= 0);
  const maxSeconds = Math.max(0, ...values);
  for (const stage of stages) {
    const item = document.querySelector(`[data-timing="${stage}"]`);
    if (!item) continue;
    item.querySelector('strong').textContent = formatQueryDuration(result, stage);
    const seconds = result[`${stage}_seconds`];
    const hasTiming = typeof seconds === 'number' && Number.isFinite(seconds) && seconds >= 0;
    item.classList.toggle('no-timing', !hasTiming);
    const ratio = hasTiming && maxSeconds > 0 ? seconds / maxSeconds : 0;
    item.style.setProperty('--timing-width', `${ratio * 100}%`);
  }
  const raster = document.querySelector('[data-timing="rasterization"]');
  if (raster) raster.classList.toggle('hidden', result.dimension !== '2d');
}
function fitResultCount() {
  const count = $('result-count');
  const text = count.textContent.trim();
  const style = getComputedStyle(count);
  const maxSize = 27;
  const minSize = 8;
  const context = document.createElement('canvas').getContext('2d');
  context.font = style.font.replace(/\d+(?:\.\d+)?px/, `${maxSize}px`);
  const measuredWidth = context.measureText(text).width;
  const availableWidth = count.clientWidth;
  if (!measuredWidth || !availableWidth) return;
  const fittedSize = Math.min(maxSize, maxSize * availableWidth / measuredWidth);
  count.style.fontSize = `${Math.max(minSize, fittedSize)}px`;
  count.style.whiteSpace = 'nowrap';
  count.style.overflowWrap = 'normal';
}
const resultCountResizeObserver = new ResizeObserver(fitResultCount);
resultCountResizeObserver.observe($('result-count').parentElement);
function showResult(result) {
  visualizationCleanup?.();
  visualizationCleanup = null;
  visualizationKey = null;
  lastResult = result;
  $('empty-state').classList.add('hidden');
  $('result-state').classList.remove('hidden');
  $('result-status').textContent = result.ok ? '完成' : `失败 (${result.exit_code})`;
  $('result-status').className = result.ok ? 'success' : 'failure';
  $('result-count').textContent = result.count == null ? '—' : result.count;
  requestAnimationFrame(fitResultCount);
  updateTiming(result);
  $('result-command').textContent = result.command_text || result.command.map(part =>
    /^[A-Za-z0-9_@%+=:,./-]+$/.test(part) ? part : "'" + part.replace(/'/g, "'\"'\"'") + "'"
  ).join(' ');
  $('visualize-button').disabled = !result.ok || !result.result_id;
  $('visualize-status').textContent = '';
  updateOutput(result, 'stdout');
  activateOutput('stdout');
}
function activateOutput(outputName) {
  document.querySelectorAll('.tab').forEach(tab => tab.classList.toggle('active', tab.dataset.output === outputName));
  const visual = outputName === 'visualization';
  setVisible('output-help', !visual);
  setVisible('result-output', !visual);
  setVisible('visualization-output', visual);
  if (!visual && lastResult) updateOutput(lastResult, outputName);
  if (visual && !visualizationKey) {
    $('visualization-output').textContent = '尚未生成可视化';
  }
}
document.querySelectorAll('.tab').forEach(tab => tab.addEventListener('click', () => activateOutput(tab.dataset.output)));

$('visualize-button').addEventListener('click', async () => {
  const idText = $('visualize-id').value.trim();
  if (!lastResult?.result_id) {
    $('visualize-status').textContent = '请先完成一次查询';
    return;
  }
  if (!/^\d+$/.test(idText)) {
    $('visualize-status').textContent = '请输入非负整数 ID';
    $('visualize-id').focus();
    return;
  }
  const button = $('visualize-button');
  button.disabled = true;
  $('visualize-status').textContent = '正在读取查询结果…';
  activateOutput('visualization');
  try {
    const response = await fetch('/api/visualize', {
      method: 'POST',
      headers: {'Content-Type': 'application/json'},
      body: JSON.stringify({result_id: lastResult.result_id, id: idText}),
    });
    const data = await response.json();
    if (!response.ok) throw new Error(data.error || '可视化失败');
    const key = `${lastResult.result_id}:${idText}`;
    if (key !== visualizationKey) {
      visualizationCleanup?.();
      visualizationCleanup = null;
      const view = await import('/visualization.js?v=20261110-visualize-7');
      if (data.dimension === '2d') view.render2D($('visualization-output'), data);
      else visualizationCleanup = await view.render3D($('visualization-output'), data);
      visualizationKey = key;
    }
    $('visualize-status').textContent = `${data.shown_neighbors} 个关联对象`;
  } catch (error) {
    visualizationCleanup?.();
    visualizationCleanup = null;
    visualizationKey = null;
    const message = lastResult.dimension === '3d'
      ? '三维渲染模块加载失败，请检查浏览器的网络连接。'
      : `可视化失败：${error.message}`;
    $('visualization-output').textContent = message;
    $('visualize-status').textContent = '读取失败';
  } finally {
    button.disabled = false;
  }
});
$('visualize-id').addEventListener('keydown', event => {
  if (event.key === 'Enter') {
    event.preventDefault();
    $('visualize-button').click();
  }
});

form.addEventListener('submit', async event => {
  event.preventDefault();
  const button = form.querySelector('.run-button');
  button.disabled = true;
  button.classList.add('running');
  $('run-label').textContent = '正在运行…';
  try {
    const dimension = selected('dimension');
    const payload = {
      dimension, backend: selected('backend'), query: query.value,
      source: $('source').value, target: $('target').value, distance: $('distance').value,
      target_object: $('target-object').value,
      k: $('k').value, threads: $('threads').value,
      progressive_filter: dimension === '3d' && $('progressive-filter').checked,
      vpr: $('vpr').value, bitwidth: $('bitwidth').value, 'compute_threads': $('compute-threads').value,
      hierarchy: dimension === '2d' && $('hierarchy').checked,
      approximation: dimension === '2d' && $('approximation').checked,
      print_result: true,
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
