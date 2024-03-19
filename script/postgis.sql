create table nuclei_box(id integer primary key, hausdorff float, phausdorff float, geom geometry(GeometryZ));
create table nuclei_20(id integer primary key, hausdorff float, phausdorff float, geom geometry(GeometryZ));
create table nuclei_40(id integer primary key, hausdorff float, phausdorff float, geom geometry(GeometryZ));
create table nuclei_60(id integer primary key, hausdorff float, phausdorff float, geom geometry(GeometryZ));
create table nuclei_80(id integer primary key, hausdorff float, phausdorff float, geom geometry(GeometryZ));
create table nuclei_100(id integer primary key, hausdorff float, phausdorff float, geom geometry(GeometryZ));

create table vessel_box(id integer primary key, hausdorff float, phausdorff float, geom geometry(GeometryZ));
create table vessel_20(id integer primary key, hausdorff float, phausdorff float, geom geometry(GeometryZ));
create table vessel_40(id integer primary key, hausdorff float, phausdorff float, geom geometry(GeometryZ));
create table vessel_60(id integer primary key, hausdorff float, phausdorff float, geom geometry(GeometryZ));
create table vessel_80(id integer primary key, hausdorff float, phausdorff float, geom geometry(GeometryZ));
create table vessel_100(id integer primary key, hausdorff float, phausdorff float, geom geometry(GeometryZ));


insert into nuclei_box select id, hausdorff, phausdorff, ST_3DExtent(geom) from nuclei_100 group by id;
insert into vessel_box select id, hausdorff, phausdorff, ST_3DExtent(geom) from vessel_100 group by id;

CREATE INDEX idx_3d_geom_gist ON my_table USING GIST (geom);        
