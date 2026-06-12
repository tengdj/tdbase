#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <limits>
#include <string>
#include <cstdlib>
#include <float.h>
#include <assert.h>
#include <queue>

// 编译推荐：-std=c++17 -O3 -mavx2 -ffast-math -flto

// ===================== 基础数学结构 =====================
struct Vec3 {
    float x, y, z;

    Vec3() : x(0), y(0), z(0) {}
    Vec3(float x_, float y_, float z_) : x(x_), y(y_), z(z_) {}

    inline Vec3 operator+(const Vec3& v) const { return {x + v.x, y + v.y, z + v.z}; }
    inline Vec3 operator-(const Vec3& v) const { return {x - v.x, y - v.y, z - v.z}; }
    inline Vec3 operator*(float s) const { return {x * s, y * s, z * s}; }
    inline Vec3 operator/(float s) const { return {x / s, y / s, z / s}; }

    inline float dot(const Vec3& v) const { return x*v.x + y*v.y + z*v.z; }
    inline Vec3 cross(const Vec3& v) const {
        return { y*v.z - z*v.y, z*v.x - x*v.z, x*v.y - y*v.x };
    }
    inline float length2() const { return x*x + y*y + z*z; }
    inline float length() const { return std::sqrt(length2()); }
};

std::ostream& operator<<(std::ostream& os, const Vec3& v)
{
    os << "(" << v.x << ", " << v.y << ", " << v.z << ")";
    return os;
}

struct AABB {
    Vec3 min, max;
    AABB() {
        min = Vec3(FLT_MAX, FLT_MAX, FLT_MAX);
        max = Vec3(-FLT_MAX, -FLT_MAX, -FLT_MAX);
    }

    void expand(const Vec3& p) {
        min.x = std::min(min.x, p.x);
        min.y = std::min(min.y, p.y);
        min.z = std::min(min.z, p.z);
        max.x = std::max(max.x, p.x);
        max.y = std::max(max.y, p.y);
        max.z = std::max(max.z, p.z);
    }

    void expand(const AABB& box) {
        expand(box.min);
        expand(box.max);
    }

    // 点到AABB 最小距离平方（内部/边界返回0）
    inline float distance2(const Vec3& p) const
    {
        float dx = 0.0f;
        if (p.x < min.x) dx = min.x - p.x;
        else if (p.x > max.x) dx = p.x - max.x;

        float dy = 0.0f;
        if (p.y < min.y) dy = min.y - p.y;
        else if (p.y > max.y) dy = p.y - max.y;

        float dz = 0.0f;
        if (p.z < min.z) dz = min.z - p.z;
        else if (p.z > max.z) dz = p.z - max.z;

        return dx*dx + dy*dy + dz*dz;
    }

    // 判断点是否在AABB内(含边界)
    inline bool contains(const Vec3& p) const
    {
        return (p.x >= min.x && p.x <= max.x) &&
               (p.y >= min.y && p.y <= max.y) &&
               (p.z >= min.z && p.z <= max.z);
    }

    // 点到AABB 最大距离平方（点在内部时使用）
    inline float maxDistance2(const Vec3& p) const
    {
		float dx = std::max(fabs(p.x - min.x), fabs(p.x - max.x));
		float dy = std::max(fabs(p.y - min.y), fabs(p.y - max.y));
		float dz = std::max(fabs(p.z - min.z), fabs(p.z - max.z));
		return dx*dx + dy*dy + dz*dz;
    }

    // 将当前AABB包围盒导出为立方体OFF网格
	bool exportToOFF(const std::string& path) const
	{
		std::ofstream f(path);
		if (!f.is_open())
			return false;

		// AABB 8个顶点
		Vec3 v[8] = {
			{min.x, min.y, min.z},
			{max.x, min.y, min.z},
			{max.x, max.y, min.z},
			{min.x, max.y, min.z},
			{min.x, min.y, max.z},
			{max.x, min.y, max.z},
			{max.x, max.y, max.z},
			{min.x, max.y, max.z}
		};

		// 6个矩形面，每个拆2个三角，共12个三角面
		int faces[12][3] = {
			{0,1,2}, {0,2,3}, // 底面 z-min
			{4,6,5}, {4,7,6}, // 顶面 z-max
			{0,3,7}, {0,7,4}, // 左面 x-min
			{1,5,6}, {1,6,2}, // 右面 x-max
			{0,4,5}, {0,5,1}, // 前面 y-min
			{3,2,6}, {3,6,7}  // 后面 y-max
		};

		f << "OFF\n";
		f << "8 12 0\n";

		// 写入8个顶点
		for (int i = 0; i < 8; ++i)
		{
			f << v[i].x << " " << v[i].y << " " << v[i].z << "\n";
		}

		// 写入12个三角面片
		for (int i = 0; i < 12; ++i)
		{
			f << "3 " << faces[i][0] << " " << faces[i][1] << " " << faces[i][2] << "\n";
		}

		f.close();
		return true;
	}

};

// 输出 AABB: [min, max]
std::ostream& operator<<(std::ostream& os, const AABB& box)
{
    os << "AABB{ min:" << box.min << ", max:" << box.max << " }";
    return os;
}

struct Triangle {
    int idx;
    Vec3 v0, v1, v2;
    AABB aabb;
    Vec3 centroid;

    Triangle(int id, const Vec3& a, const Vec3& b, const Vec3& c)
        : idx(id), v0(a), v1(b), v2(c)
    {
        aabb.expand(v0);
        aabb.expand(v1);
        aabb.expand(v2);
        centroid = (v0 + v1 + v2) * (1.0f / 3.0f);
    }
    bool exportToOFF(const std::string& path) const
	{
		std::ofstream f(path);
		if (!f.is_open())
			return false;

		f << "OFF\n";
		// 3顶点 1面片 0边
		f << "3 1 0\n";
		// 输出三个顶点
		f << v0.x << " " << v0.y << " " << v0.z << "\n";
		f << v1.x << " " << v1.y << " " << v1.z << "\n";
		f << v2.x << " " << v2.y << " " << v2.z << "\n";
		// 三角面：顶点数 + 三个索引
		f << "3 0 1 2\n";

		f.close();
		return true;
	}
};

// 输出 Triangle 信息
std::ostream& operator<<(std::ostream& os, const Triangle& tri)
{
    os << "Triangle{ idx:" << tri.idx
       << ", v0:" << tri.v0
       << ", v1:" << tri.v1
       << ", v2:" << tri.v2
       << ", aabb:" << tri.aabb
       << " }";
    return os;
}

struct BVHNode {
    AABB aabb;
    int left, right;
    int triStart, triCnt;
    BVHNode() : left(-1), right(-1), triStart(0), triCnt(0) {}
};

// BVH遍历单元：节点下标 + 该节点到查询点的最小距离平方
struct BVHTraversalNode
{
    int nodeIdx;
    float minDist2;

    BVHTraversalNode(int idx, float d2) : nodeIdx(idx), minDist2(d2) {}

    // 小顶堆：最小距离优先
    bool operator>(const BVHTraversalNode& other) const
    {
        return minDist2 > other.minDist2;
    }
};

struct PointQueryRes {
    float dist;
    int nearTriID;
    PointQueryRes() : dist(FLT_MAX), nearTriID(-1) {}
};

std::ostream& operator<<(std::ostream& os, const PointQueryRes& res)
{
	os << "PointQueryRes{ dist: " << res.dist
	   << ", nearTriID: " << res.nearTriID << " }";
	return os;
}

// 面片统计结果
struct TriStatResult {
    // A面片自身采样点 -> B网格 距离最大值
    float maxSelfSampleToOther;
    // 归属该A面片的所有B采样点 -> A面片 距离最大值
    float maxOtherSampleToSelf;

    TriStatResult()
        : maxSelfSampleToOther(0.0f),
          maxOtherSampleToSelf(0.0f)
    {}
};

// ===================== 网格类 =====================
struct Mesh {
    std::vector<Vec3> vertices;
    std::vector<Triangle> triangles;
    std::vector<BVHNode> bvh;
    std::vector<int> triIndices;

    bool loadOFF(const std::string& filePath);
    void buildBVH();
    // 【修改】按固定空间间隔采样，替代按点数采样
    std::vector<Vec3> sampleByStep(float step) const;
    PointQueryRes queryNearest(const Vec3& p) const;
    PointQueryRes queryNearestBruteForce(const Vec3& p) const;
};

// ===================== 点-三角形距离 =====================
inline float pointTriangleDistance2(const Vec3& p, const Triangle& tri)
{
    const Vec3& a = tri.v0;
    const Vec3& b = tri.v1;
    const Vec3& c = tri.v2;

    Vec3 ab = b - a;
    Vec3 ac = c - a;
    Vec3 ap = p - a;

    float d1 = ab.dot(ap);
    float d2 = ac.dot(ap);
    if (d1 <= 0.0f && d2 <= 0.0f)
        return ap.length2();

    Vec3 bp = p - b;
    float d3 = ab.dot(bp);
    float d4 = ac.dot(bp);
    if (d3 >= 0.0f && d4 <= d3)
        return bp.length2();

    float vc = d1 * d4 - d3 * d2;
    if (vc <= 0.0f && d1 >= 0.0f && d3 <= 0.0f)
    {
        float v = d1 / (d1 - d3);
        return (ap + ab * (-v)).length2();
    }

    Vec3 cp = p - c;
    float d5 = ab.dot(cp);
    float d6 = ac.dot(cp);
    if (d6 >= 0.0f && d5 <= d6)
        return cp.length2();

    float vb = d5 * d2 - d1 * d6;
    if (vb <= 0.0f && d2 >= 0.0f && d6 <= 0.0f)
    {
        float w = d2 / (d2 - d6);
        return (ap + ac * (-w)).length2();
    }

    float va = d3 * d6 - d5 * d4;
    if (va <= 0.0f && (d3 - d5) >= 0.0f && (d4 - d6) >= 0.0f)
    {
        float u = (d3 - d5) / ((d3 - d5) + (d4 - d6));
        return (bp + (c - b) * (-u)).length2();
    }

    float denom = 1.0f / (va + vb + vc);
    float v = vb * denom;
    float w = vc * denom;
    return (ap + ab * (-v) + ac * (-w)).length2();
}

inline float pointTriangleDistance(const Vec3& p, const Triangle& tri)
{
    return std::sqrt(pointTriangleDistance2(p, tri));
}


// ===================== BVH 构建与查询 =====================
void buildBVHRecursive(std::vector<BVHNode>& nodes,
                       std::vector<int>& triIndices,
                       int nodeIdx, int start, int end,
                       const std::vector<Triangle>& tris)
{
    BVHNode& node = nodes[nodeIdx];
    node.triStart = start;
    node.triCnt = end - start;
    node.aabb = AABB();

    for (int i = start; i < end; ++i){
    	node.aabb.expand(tris[triIndices[i]].aabb);
    }

    if (node.triCnt < 4)
    {
        node.left = node.right = -1;
        return;
    }

    Vec3 diag = node.aabb.max - node.aabb.min;
    int axis = 0;
    if (diag.y > diag.x) axis = 1;
    if (diag.z > (axis == 0 ? diag.x : diag.y)) axis = 2;

    float midVal = ((axis == 0 ? node.aabb.min.x : (axis == 1 ? node.aabb.min.y : node.aabb.min.z))
                  + (axis == 0 ? node.aabb.max.x : (axis == 1 ? node.aabb.max.y : node.aabb.max.z))) * 0.5f;

    auto midIt = std::partition(triIndices.begin() + start, triIndices.begin() + end,
        [&](int t) {
            const Vec3& c = tris[t].centroid;
            float val = (axis == 0) ? c.x : (axis == 1 ? c.y : c.z);
            return val < midVal;
        });

    int midIdx = midIt - triIndices.begin();
    if (midIdx == start || midIdx == end)
        midIdx = (start + end) / 2;

    node.left = (int)nodes.size();
    nodes.emplace_back();
    node.right = (int)nodes.size();
    nodes.emplace_back();

    buildBVHRecursive(nodes, triIndices, node.left, start, midIdx, tris);
    buildBVHRecursive(nodes, triIndices, node.right, midIdx, end, tris);
}

void bvhQueryOpt(
    const std::vector<BVHNode>& nodes,
    const std::vector<Triangle>& tris,
    const std::vector<int> triIndices,
    const Vec3& p,
    float& globalMinDist2,
    int& nearestTriId)
{
    // 优先队列：小顶堆，离点越近的节点优先处理
    std::priority_queue<BVHTraversalNode,
                        std::vector<BVHTraversalNode>,
                        std::greater<BVHTraversalNode>> pq;

    // 根节点入队
    float rootMinD2 = nodes[0].aabb.distance2(p);
    pq.emplace(0, rootMinD2);

    while (!pq.empty())
    {
        auto curr = pq.top();
        pq.pop();

        int nodeIdx = curr.nodeIdx;
        float nodeMinD2 = curr.minDist2;
        const BVHNode& node = nodes[nodeIdx];

        printf("%d\t%f %f\n",nodeIdx,nodeMinD2,globalMinDist2);

        // ========== 剪枝：当前节点最小可能距离 > 全局已知最近距离，直接跳过 ==========
        if (nodeMinD2 > globalMinDist2)
            continue;
        // 叶子节点：遍历内部所有面片
        if (node.left == -1)
        {
        	//printf("%d %d %d \n",nodeIdx,node.triStart,node.triCnt);
            for (int i = node.triStart; i < node.triStart + node.triCnt; ++i)
            {
                const Triangle& tri = tris[triIndices[i]];
                float d2 = pointTriangleDistance2(p, tri);
                //printf("%f in [%f %f]\n",d2,node.aabb.distance2(p),node.aabb.maxDistance2(p));
                if (d2 < globalMinDist2)
                {
                    globalMinDist2 = d2;
                    nearestTriId = tri.idx;
                	std::cout<<"get here"<<std::sqrt(d2)<<std::endl;
                }
                //std::cout<<tri<<std::endl;
//                char path[256];
//                sprintf(path,"tri%d.off",i);
//                tri.exportToOFF(path);
            }
//            std::cout<<node.aabb<<std::endl;
//            node.aabb.exportToOFF("aabb.off");
//            exit(0);
            continue;
        }

        // 非叶子节点：计算左右子节点距离区间，入队
        const BVHNode& leftNode  = nodes[node.left];
        const BVHNode& rightNode = nodes[node.right];

        globalMinDist2 = std::min(globalMinDist2,leftNode.aabb.maxDistance2(p));
        globalMinDist2 = std::min(globalMinDist2,rightNode.aabb.maxDistance2(p));

        // 左子节点 最小距离平方
        float leftMinD2 = leftNode.aabb.distance2(p);
        // 右子节点 最小距离平方
        float rightMinD2 = rightNode.aabb.distance2(p);

        // 仅当子节点有机会更优时才入队（前置剪枝）
        if (leftMinD2 <= globalMinDist2)
            pq.emplace(node.left, leftMinD2);
        if (rightMinD2 <= globalMinDist2)
            pq.emplace(node.right, rightMinD2);
    }
}

// 输入：BVH 节点数组，输出路径；将所有节点 AABB 导出为 OFF 网格
bool exportBVHAABBToOFF(const std::vector<BVHNode>& bvhNodes, const std::string& filePath)
{
    if (bvhNodes.empty())
        return false;

    std::vector<Vec3> vertices;
    std::vector<std::vector<int>> faces;

    // 单个 AABB 生成 8 顶点 + 12 三角面片
    auto addAABBGeometry = [&](const AABB& box)
    {
        int baseIdx = static_cast<int>(vertices.size());
        // AABB 八个顶点
        vertices.emplace_back(box.min.x, box.min.y, box.min.z);
        vertices.emplace_back(box.max.x, box.min.y, box.min.z);
        vertices.emplace_back(box.max.x, box.max.y, box.min.z);
        vertices.emplace_back(box.min.x, box.max.y, box.min.z);

        vertices.emplace_back(box.min.x, box.min.y, box.max.z);
        vertices.emplace_back(box.max.x, box.min.y, box.max.z);
        vertices.emplace_back(box.max.x, box.max.y, box.max.z);
        vertices.emplace_back(box.min.x, box.max.y, box.max.z);

        // 6个面，每个面拆2个三角形
        // 底面 z=min
        faces.push_back({3, baseIdx+0, baseIdx+1, baseIdx+2});
        faces.push_back({3, baseIdx+0, baseIdx+2, baseIdx+3});
        // 顶面 z=max
        faces.push_back({3, baseIdx+4, baseIdx+6, baseIdx+5});
        faces.push_back({3, baseIdx+4, baseIdx+7, baseIdx+6});
        // 左面 x=min
        faces.push_back({3, baseIdx+0, baseIdx+3, baseIdx+7});
        faces.push_back({3, baseIdx+0, baseIdx+7, baseIdx+4});
        // 右面 x=max
        faces.push_back({3, baseIdx+1, baseIdx+5, baseIdx+6});
        faces.push_back({3, baseIdx+1, baseIdx+6, baseIdx+2});
        // 前面 y=min
        faces.push_back({3, baseIdx+0, baseIdx+4, baseIdx+5});
        faces.push_back({3, baseIdx+0, baseIdx+5, baseIdx+1});
        // 后面 y=max
        faces.push_back({3, baseIdx+3, baseIdx+2, baseIdx+6});
        faces.push_back({3, baseIdx+3, baseIdx+6, baseIdx+7});
    };

    // 遍历所有 BVH 节点
    for (const auto& node : bvhNodes)
    {
        addAABBGeometry(node.aabb);
    }

    // 写入标准 OFF 文件
    std::ofstream fout(filePath);
    if (!fout.is_open())
        return false;

    fout << "OFF\n";
    fout << vertices.size() << " " << faces.size() << " 0\n";

    // 写入顶点
    for (const auto& v : vertices)
    {
        fout << v.x << " " << v.y << " " << v.z << "\n";
    }
    // 写入三角面片
    for (const auto& f : faces)
    {
        fout << f[0] << " " << f[1] << " " << f[2] << " " << f[3] << "\n";
    }

    fout.close();
    return true;
}

int bp = 0;
void Mesh::buildBVH()
{
    bvh.clear();
    bvh.reserve(triangles.size() * 2);
    bvh.emplace_back();
    triIndices.resize(triangles.size());
    for (int i = 0; i < (int)triangles.size(); ++i)
    	triIndices[i] = i;

    buildBVHRecursive(bvh, triIndices, 0, 0, (int)triangles.size(), triangles);
//    char path[256];
//    sprintf(path,"../../data/bvh%d.OFF",bp++);
//    exportBVHAABBToOFF(bvh, path);
}

PointQueryRes Mesh::queryNearest(const Vec3& p) const
{
    PointQueryRes res;
    float minD2 = FLT_MAX;
    int nearTri = -1;
    bvhQueryOpt(bvh, triangles, triIndices, p, minD2, nearTri);

    res.dist = std::sqrt(minD2);
    res.nearTriID = nearTri;
    return res;
}

PointQueryRes  Mesh::queryNearestBruteForce(const Vec3& p) const
{
    PointQueryRes res;
    float minDist2 = FLT_MAX;
    int nearTriId = -1;

    for (const auto& tri : triangles)
    {
        float d2 = pointTriangleDistance2(p, tri);
        if (d2 < minDist2)
        {
            minDist2 = d2;
            nearTriId = tri.idx;
        }
    }

    res.dist = std::sqrt(minDist2);
    res.nearTriID = nearTriId;

//    printf("%d\t%f\n",res.nearTriID, res.dist);
//    exit(0);
    return res;
}

// ===================== OFF 读取 =====================
bool Mesh::loadOFF(const std::string& filePath)
{
    std::ifstream fin(filePath);
    if (!fin.is_open())
    {
        std::cerr << "Open file failed: " << filePath << std::endl;
        return false;
    }

    std::string line;
    std::getline(fin, line);
    while (line.empty()) std::getline(fin, line);

    int vCnt, fCnt, eCnt;
    fin >> vCnt >> fCnt >> eCnt;

    vertices.reserve(vCnt);
    for (int i = 0; i < vCnt; ++i)
    {
        float x, y, z;
        fin >> x >> y >> z;
        vertices.emplace_back(x, y, z);
    }

    triangles.reserve(fCnt);
    for (int i = 0; i < fCnt; ++i)
    {
        int n, v0, v1, v2;
        fin >> n >> v0 >> v1 >> v2;
        if (n == 3)
        {
            triangles.emplace_back(i, vertices[v0], vertices[v1], vertices[v2]);
        }
    }
    fin.close();
    buildBVH();
    return true;
}

// ===================== 【核心修改】按固定空间间隔采样 =====================
// 在三角形内部，按指定空间步长 step 生成等距规则采样点
std::vector<Vec3> Mesh::sampleByStep(float step) const
{
    std::vector<Vec3> points;
    if (step <= 1e-6f || triangles.empty())
        return points;

    // 先预估容量（粗略）
    for (const auto& tri : triangles)
    {
        Vec3 e1 = tri.v1 - tri.v0;
        Vec3 e2 = tri.v2 - tri.v0;
        float len1 = e1.length();
        float len2 = e2.length();
        int cnt1 = static_cast<int>(len1 / step) + 2;
        int cnt2 = static_cast<int>(len2 / step) + 2;
        points.reserve(points.capacity() + cnt1 * cnt2 / 2);
    }

    // 遍历每个三角面片采样
    for (const auto& tri : triangles)
    {
        const Vec3& v0 = tri.v0;
        const Vec3& v1 = tri.v1;
        const Vec3& v2 = tri.v2;

        Vec3 edgeU = v1 - v0;
        Vec3 edgeV = v2 - v0;

        float lenU = edgeU.length();
        float lenV = edgeV.length();

        // 计算两条边上的分段数
        int numU = static_cast<int>(lenU / step) + 1;
        int numV = static_cast<int>(lenV / step) + 1;

        numU = std::max(numU, 1);
        numV = std::max(numV, 1);

        // 重心坐标均匀步进，保证空间距离近似等于 step
        for (int i = 0; i <= numU; ++i)
        {
            float u = static_cast<float>(i) / numU;
            for (int j = 0; j <= numV; ++j)
            {
                float v = static_cast<float>(j) / numV;
                // 仅保留三角形内部点：u + v <= 1
                if (u + v > 1.0f)
                    continue;

                float w = 1.0f - u - v;
                Vec3 p = v0 * w + v1 * u + v2 * v;
                points.push_back(p);
            }
        }
    }
    return points;
}

// 将点云向量导出为 OFF 点云文件
bool exportPointsToOFF(const std::vector<Vec3>& points, const std::string& outPath)
{
    std::ofstream fout(outPath);
    if (!fout.is_open())
    {
        std::cerr << "Failed to open output file: " << outPath << std::endl;
        return false;
    }

    // OFF 头部: 顶点数 面片数 边数
    fout << "OFF\n";
    fout << points.size() << " 0 0\n";

    // 逐行写入坐标
    for (const auto& p : points)
    {
        fout << p.x << " " << p.y << " " << p.z << "\n";
    }

    fout.close();
    return true;
}
// ===================== 主业务接口 =====================
/**
 * @param pathA 网格A路径
 * @param pathB 网格B路径
 * @param sampleStep 全局固定采样间隔（空间距离）
 * @return A每个面片的两类距离统计值
 */
std::vector<TriStatResult> computeMeshStat(
    const std::string& pathA,
    const std::string& pathB,
    float sampleStep)
{
    Mesh meshA, meshB;
    if (!meshA.loadOFF(pathA) || !meshB.loadOFF(pathB))
        return {};

    const int triNumA = (int)meshA.triangles.size();
    std::vector<TriStatResult> result(triNumA);

    // ========== 1. A网格按固定间隔采样，并记录每个采样点所属A面片 ==========
    std::vector<Vec3> sampleA;
    std::vector<int> sample2TriA;

    // 逐面片采样 + 归属标记
    for (int triId = 0; triId < triNumA; ++triId)
    {
        const Triangle& tri = meshA.triangles[triId];
        const Vec3& v0 = tri.v0;
        const Vec3& v1 = tri.v1;
        const Vec3& v2 = tri.v2;

        Vec3 edgeU = v1 - v0;
        Vec3 edgeV = v2 - v0;
        float lenU = edgeU.length();
        float lenV = edgeV.length();

        int numU = static_cast<int>(lenU / sampleStep) + 1;
        int numV = static_cast<int>(lenV / sampleStep) + 1;
        numU = std::max(numU, 1);
        numV = std::max(numV, 1);

        for (int i = 0; i <= numU; ++i)
        {
            float u = static_cast<float>(i) / numU;
            for (int j = 0; j <= numV; ++j)
            {
                float v = static_cast<float>(j) / numV;
                if (u + v > 1.0f) continue;
                float w = 1.0f - u - v;
                Vec3 p = v0 * w + v1 * u + v2 * v;
                sampleA.push_back(p);
                sample2TriA.push_back(triId);
            }
        }
    }

    // 统计：A各面片采样点到B的距离最大值
    for (size_t i = 0; i < sampleA.size(); ++i)
    {
        PointQueryRes res = meshB.queryNearest(sampleA[i]);
        PointQueryRes res1 = meshB.queryNearestBruteForce(sampleA[i]);
        std::cout<<res<<std::endl;
        std::cout<<res1<<std::endl;
        exit(0);
        int triId = sample2TriA[i];
        if (res.dist > result[triId].maxSelfSampleToOther)
        {
            result[triId].maxSelfSampleToOther = res.dist;
        }
    }

    // ========== 2. B网格按相同间隔采样，归类到A最近面片并统计最大值 ==========
    std::vector<Vec3> sampleB = meshB.sampleByStep(sampleStep);
    for (const Vec3& p : sampleB)
    {
        //PointQueryRes res = meshA.queryNearest(p);
        PointQueryRes res = meshA.queryNearestBruteForce(p);
        int belongTriA = res.nearTriID;
        assert(belongTriA >= 0 && belongTriA < triNumA);
        if (res.dist > result[belongTriA].maxOtherSampleToSelf)
        {
            result[belongTriA].maxOtherSampleToSelf = res.dist;
        }
    }

    exportPointsToOFF(sampleB,"point.off");

    std::cout<<sampleStep<<" "<<sampleA.size()<<" "<<sampleB.size()<<std::endl;
    return result;
}

// ===================== 测试入口 =====================
int main(int argc, char **argv)
{
    std::string fileA = argv[1];
    std::string fileB = argv[2];
    float step = atof(argv[3]);

    auto stats = computeMeshStat(fileA, fileB, step);

    // 输出结果
    for (int i = 0; i < (int)stats.size(); ++i)
    {
        std::cout << "Tri " << i
                  << " | Max(A_sample->B): " << stats[i].maxSelfSampleToOther
                  << " | Max(B_sample->A_tri): " << stats[i].maxOtherSampleToSelf
                  << std::endl;
    }

    return 0;
}