const SVG_NS = 'http://www.w3.org/2000/svg';
const TARGET_COLOR = '#f05252';
const BOUNDARY_COLOR = '#f1c45a';
function colorFor(index) {
  const hue = ((index * 137.508 + 170) % 360) / 60;
  const saturation = 0.66;
  const lightness = 0.53;
  const chroma = (1 - Math.abs(2 * lightness - 1)) * saturation;
  const secondary = chroma * (1 - Math.abs(hue % 2 - 1));
  const base = lightness - chroma / 2;
  const channels = hue < 1 ? [chroma, secondary, 0]
    : hue < 2 ? [secondary, chroma, 0]
      : hue < 3 ? [0, chroma, secondary]
        : hue < 4 ? [0, secondary, chroma]
          : hue < 5 ? [secondary, 0, chroma] : [chroma, 0, secondary];
  return `#${channels.map(value => Math.round((value + base) * 255).toString(16).padStart(2, '0')).join('')}`;
}

function element(name, attributes = {}) {
  const node = document.createElementNS(SVG_NS, name);
  for (const [key, value] of Object.entries(attributes)) node.setAttribute(key, String(value));
  return node;
}

function mapCenter(item) {
  return item.center;
}

function polygonCentroid(ring) {
  let twiceArea = 0;
  let x = 0;
  let y = 0;
  for (let i = 0; i < ring.length; i += 1) {
    const [ax, ay] = ring[i];
    const [bx, by] = ring[(i + 1) % ring.length];
    const cross = ax * by - bx * ay;
    twiceArea += cross;
    x += (ax + bx) * cross;
    y += (ay + by) * cross;
  }
  if (Math.abs(twiceArea) < 1e-12) {
    const xs = ring.map(point => point[0]);
    const ys = ring.map(point => point[1]);
    return [(Math.min(...xs) + Math.max(...xs)) / 2, (Math.min(...ys) + Math.max(...ys)) / 2];
  }
  return [x / (3 * twiceArea), y / (3 * twiceArea)];
}

function pathForRings(rings) {
  return rings.map(ring => ring.map((point, index) => `${index ? 'L' : 'M'}${point[0]} ${point[1]}`).join(' ') + ' Z').join(' ');
}

function pathForPoints(points) {
  return points.map((point, index) => `${index ? 'L' : 'M'}${point[0]} ${point[1]}`).join(' ') + ' Z';
}

function appendSvgIdLabel(drawing, id, center, color, size) {
  const label = String(id);
  const width = size * (label.length * 0.68 + 0.9);
  const height = size * 1.55;
  const group = element('g', {
    class: 'visual-id-label',
    transform: `translate(${center[0]} ${center[1] + height * 0.9}) scale(1 -1)`,
  });
  group.append(element('rect', {
    x: -width / 2, y: -height / 2, width, height, rx: height * 0.3,
    fill: '#ffffff', 'fill-opacity': '0.94', stroke: color, 'stroke-width': size * 0.09,
  }));
  const text = element('text', {
    x: 0, y: 0, fill: '#173238', 'font-size': size * 0.88,
    'font-family': 'system-ui, sans-serif', 'font-weight': 700,
    'text-anchor': 'middle', 'dominant-baseline': 'central',
  });
  text.textContent = label;
  group.append(text);
  drawing.append(group);
}

function addLegend(container, data) {
  const legend = document.createElement('div');
  legend.className = 'visualization-legend';
  const items = [[`Target · ${data.target.id}`, TARGET_COLOR], ...data.neighbors.map((item, index) => [`ID ${item.id}`, colorFor(index)])];
  if (data.dimension === '3d' && data.distance != null) items.push([`DWithin ${data.distance} ${data.distance_unit}`, BOUNDARY_COLOR]);
  for (const [label, color] of items) {
    const item = document.createElement('span');
    item.className = 'visualization-legend-item';
    const swatch = document.createElement('i');
    swatch.style.setProperty('--legend-color', color);
    item.append(swatch, document.createTextNode(label));
    legend.append(item);
  }
  container.append(legend);
}

function addSummary(container, data) {
  const summary = document.createElement('div');
  summary.className = 'visualization-summary';
  let text = `对象 ${data.target.id} 关联 ${data.total_neighbors.toLocaleString()} 个对象，当前显示 ${data.shown_neighbors.toLocaleString()} 个。`;
  if (data.truncated) text += ' 关联对象较多，图中仅展示前 150 个。';
  if (data.simplified) text += ' 多边形轮廓已抽稀显示。';
  if (data.dimension === '3d' && data.distance != null) text += ' 超出距离球的结果点沿连线投影到球面上。';
  summary.textContent = text;
  container.append(summary);
}

export function projectWithinResultPoint(targetCenter, resultCenter, distance) {
  if (!(distance > 0)) return resultCenter;
  const offset = resultCenter.map((value, axis) => value - targetCenter[axis]);
  const centerDistance = Math.hypot(...offset);
  if (centerDistance <= distance || centerDistance === 0) return resultCenter;
  const scale = distance / centerDistance;
  return targetCenter.map((value, axis) => value + offset[axis] * scale);
}

export function render2D(container, data) {
  container.replaceChildren();
  const figure = document.createElement('div');
  figure.className = 'visualization-figure visualization-2d';
  const svg = element('svg', {role: 'img', 'aria-label': `二维对象 ${data.target.id} 的查询关联图`, preserveAspectRatio: 'xMidYMid meet'});
  const targetCenter = polygonCentroid(data.target.rings[0]);
  const allPoints = [
    ...data.target.rings.flat(),
    ...data.neighbors.flatMap(item => item.rings.flat()),
  ];
  let minX = Math.min(...allPoints.map(point => point[0]));
  let maxX = Math.max(...allPoints.map(point => point[0]));
  let minY = Math.min(...allPoints.map(point => point[1]));
  let maxY = Math.max(...allPoints.map(point => point[1]));
  const span = Math.max(maxX - minX, maxY - minY, 1e-6);
  const padding = span * 0.12;
  minX -= padding;
  maxX += padding;
  minY -= padding;
  maxY += padding;
  svg.setAttribute('viewBox', `${minX} ${-maxY} ${maxX - minX} ${maxY - minY}`);
  const drawing = element('g', {transform: 'scale(1 -1)'});

  data.neighbors.forEach((item, index) => {
    const center = polygonCentroid(item.rings[0]);
    drawing.append(element('line', {
      x1: targetCenter[0], y1: targetCenter[1], x2: center[0], y2: center[1],
      class: 'visual-link', stroke: colorFor(index),
    }));
  });
  const polygons = [data.target, ...data.neighbors];
  polygons.forEach((item, index) => {
    const target = index === 0;
    const shape = element('path', {
      d: pathForRings(item.rings),
      class: target ? 'visual-polygon visual-target' : 'visual-polygon',
      fill: target ? TARGET_COLOR : colorFor(index - 1),
      'fill-rule': 'evenodd',
      stroke: target ? '#b72e39' : colorFor(index - 1),
    });
    const title = element('title');
    title.textContent = target ? `Target ${item.id}` : `ID ${item.id}`;
    shape.append(title);
    drawing.append(shape);
  });
  const labelSize = span * 0.018;
  polygons.forEach((item, index) => {
    const center = polygonCentroid(item.rings[0]);
    appendSvgIdLabel(drawing, item.id, center, index === 0 ? TARGET_COLOR : colorFor(index - 1), labelSize);
  });
  svg.append(drawing);
  figure.append(svg);
  container.append(figure);
  addLegend(container, data);
  addSummary(container, data);
}

export async function render3D(container, data) {
  const CONTROLS_URL = '/vendor/OrbitControls.js';
  const [THREE, {OrbitControls}] = await Promise.all([import('three'), import(CONTROLS_URL)]);
  container.replaceChildren();
  const figure = document.createElement('div');
  figure.className = 'visualization-figure visualization-3d';
  const canvasHost = document.createElement('div');
  canvasHost.className = 'visualization-canvas';
  figure.append(canvasHost);
  container.append(figure);

  const renderer = new THREE.WebGLRenderer({antialias: true, alpha: false, preserveDrawingBuffer: true});
  renderer.setPixelRatio(Math.min(window.devicePixelRatio || 1, 2));
  renderer.setClearColor(0x142a30, 1);
  renderer.outputColorSpace = THREE.SRGBColorSpace;
  canvasHost.append(renderer.domElement);

  const scene = new THREE.Scene();
  scene.background = new THREE.Color('#142a30');
  const camera = new THREE.PerspectiveCamera(45, 1, 0.01, 1000);
  camera.position.set(11, 9, 13);
  scene.add(new THREE.HemisphereLight(0xffffff, 0x45616a, 2.1));
  const keyLight = new THREE.DirectionalLight(0xffffff, 2.4);
  keyLight.position.set(8, 12, 9);
  scene.add(keyLight);
  const controls = new OrbitControls(camera, renderer.domElement);
  controls.enableDamping = true;
  controls.target.set(0, 0, 0);

  const targetCenter = data.target.center;
  const distance = data.distance || 0;
  const bounds = data.target.bounds;
  const targetRadius = Math.hypot(...bounds[1].map((value, axis) => (value - bounds[0][axis]) / 2));
  const neighborOffsets = data.neighbors.map(item => {
    const center = data.distance != null
      ? projectWithinResultPoint(targetCenter, item.center, distance)
      : item.center;
    return center.map((value, axis) => value - targetCenter[axis]);
  });
  const maxNeighborDistance = Math.max(...neighborOffsets.map(offset => Math.hypot(...offset)), 0);
  const fitExtent = data.distance != null && distance > 0
    ? distance
    : Math.max(maxNeighborDistance, targetRadius, 1e-6);
  const scale = 7 / fitExtent;
  const radius = Math.max(0.10, Math.min(0.24, fitExtent * scale * 0.027));
  const neighborPoints = neighborOffsets.map(offset => new THREE.Vector3(...offset.map(value => value * scale)));
  const targetPoint = new THREE.Vector3(0, 0, 0);
  const nodeMaterial = (color, emissive) => new THREE.MeshStandardMaterial({color, emissive, emissiveIntensity: 0.16, roughness: 0.34, metalness: 0.12});

  if (data.distance != null && distance > 0) {
    const sphere = new THREE.SphereGeometry(distance * scale, 16, 8);
    const surface = new THREE.Mesh(sphere, new THREE.MeshBasicMaterial({color: BOUNDARY_COLOR, transparent: true, opacity: 0.045, side: THREE.DoubleSide, depthWrite: false}));
    const wire = new THREE.Mesh(sphere, new THREE.MeshBasicMaterial({color: BOUNDARY_COLOR, transparent: true, opacity: 0.3, wireframe: true}));
    scene.add(surface, wire);
  }
  const targetBall = new THREE.Mesh(new THREE.SphereGeometry(radius, 28, 18), nodeMaterial(TARGET_COLOR, '#7b292c'));
  targetBall.position.copy(targetPoint);
  scene.add(targetBall);
  const addIdLabel = (id, color, position) => {
    const text = String(id);
    const canvas = document.createElement('canvas');
    canvas.width = 256;
    canvas.height = 96;
    const context = canvas.getContext('2d');
    const fontSize = Math.min(48, 208 / Math.max(1, text.length * 0.68));
    context.font = `700 ${fontSize}px system-ui, sans-serif`;
    context.textAlign = 'center';
    context.textBaseline = 'middle';
    context.fillStyle = 'rgba(255,255,255,0.96)';
    context.strokeStyle = color;
    context.lineWidth = 5;
    context.beginPath();
    context.roundRect(3, 3, 250, 90, 18);
    context.fill();
    context.stroke();
    context.fillStyle = '#173238';
    context.fillText(text, 128, 48);
    const texture = new THREE.CanvasTexture(canvas);
    texture.colorSpace = THREE.SRGBColorSpace;
    const sprite = new THREE.Sprite(new THREE.SpriteMaterial({map: texture, transparent: true, depthTest: false}));
    sprite.position.copy(position);
    sprite.position.y += radius * 3.1;
    sprite.scale.set(radius * 6.4, radius * 2.7, 1);
    scene.add(sprite);
  };
  addIdLabel(data.target.id, TARGET_COLOR, targetPoint);
  data.neighbors.forEach((item, index) => {
    const point = neighborPoints[index];
    const color = colorFor(index);
    const ball = new THREE.Mesh(new THREE.SphereGeometry(radius * 0.8, 20, 14), nodeMaterial(color, color));
    ball.position.copy(point);
    scene.add(ball);
    addIdLabel(item.id, color, point);
    const direction = point.clone().normalize();
    const lineStart = direction.clone().multiplyScalar(radius);
    const lineEnd = point.clone().addScaledVector(direction, -radius * 0.8);
    if (lineStart.distanceTo(lineEnd) > 0) {
      const geometry = new THREE.BufferGeometry().setFromPoints([lineStart, lineEnd]);
      scene.add(new THREE.Line(geometry, new THREE.LineBasicMaterial({color, transparent: true, opacity: 0.76})));
    }
  });

  const maxNodeDistance = Math.max(...neighborPoints.map(point => point.length()), 0);
  const maxRadius = Math.max(distance * scale, maxNodeDistance, targetRadius * scale) + radius * 6.6;
  const cameraDirection = new THREE.Vector3(1.55, 1.15, 1.75).normalize();
  controls.update();

  let frame = 0;
  const resize = () => {
    const width = Math.max(canvasHost.clientWidth, 1);
    const height = Math.max(canvasHost.clientHeight, 1);
    renderer.setSize(width, height, false);
    camera.aspect = width / height;
    const verticalHalfFov = THREE.MathUtils.degToRad(camera.fov / 2);
    const horizontalHalfFov = Math.atan(Math.tan(verticalHalfFov) * camera.aspect);
    const fitDistance = maxRadius / Math.sin(Math.min(verticalHalfFov, horizontalHalfFov)) * 1.1;
    camera.position.copy(cameraDirection).multiplyScalar(fitDistance);
    camera.near = Math.max(0.001, maxRadius / 1000);
    camera.far = Math.max(100, fitDistance + maxRadius * 3);
    camera.updateProjectionMatrix();
    controls.update();
  };
  const observer = new ResizeObserver(resize);
  observer.observe(canvasHost);
  resize();
  const animate = () => {
    frame = requestAnimationFrame(animate);
    controls.update();
    renderer.render(scene, camera);
  };
  animate();

  addLegend(container, data);
  addSummary(container, data);
  return () => {
    cancelAnimationFrame(frame);
    observer.disconnect();
    controls.dispose();
    scene.traverse(object => {
      object.geometry?.dispose();
      object.material?.map?.dispose();
      if (Array.isArray(object.material)) object.material.forEach(material => material.dispose());
      else object.material?.dispose();
    });
    renderer.dispose();
    renderer.domElement.remove();
  };
}
