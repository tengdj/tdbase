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

CREATE INDEX idx_3d_geom_nuclei_box ON nuclei_box USING GIST (geom);        
CREATE INDEX idx_3d_geom_nuclei_20 ON nuclei_20 USING GIST (geom);        
CREATE INDEX idx_3d_geom_nuclei_40 ON nuclei_40 USING GIST (geom);        
CREATE INDEX idx_3d_geom_nuclei_60 ON nuclei_60 USING GIST (geom);        
CREATE INDEX idx_3d_geom_nuclei_80 ON nuclei_80 USING GIST (geom);        
CREATE INDEX idx_3d_geom_nuclei_100 ON nuclei_100 USING GIST (geom);        

CREATE INDEX idx_3d_geom_vessel_box ON vessel_box USING GIST (geom);        
CREATE INDEX idx_3d_geom_vessel_20 ON vessel_20 USING GIST (geom);        
CREATE INDEX idx_3d_geom_vessel_40 ON vessel_40 USING GIST (geom);        
CREATE INDEX idx_3d_geom_vessel_60 ON vessel_60 USING GIST (geom);        
CREATE INDEX idx_3d_geom_vessel_80 ON vessel_80 USING GIST (geom);        
CREATE INDEX idx_3d_geom_vessel_100 ON vessel_100 USING GIST (geom);        

insert into nuclei_box select id, hausdorff, phausdorff, ST_3DExtent(geom) from nuclei_100 group by id;
insert into vessel_box select id, hausdorff, phausdorff, ST_3DExtent(geom) from vessel_100 group by id;

