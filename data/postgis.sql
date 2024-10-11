explain analyze 

CREATE TABLE ntest(id SERIAL PRIMARY KEY, geom geometry);

select n1.id,n2.id,st_3ddistance(n1.geom, n2.geom) as dist from nuclei_3d n1, nuclei_3d n2 where n1.id<>n2.id and ST_Expand(n1.geom, 400) &&& n2.geom group n1.id order by dist asc limit 1;

select n1.id,n2.id from nuclei1 n1, nuclei2 n2 where n1.geom &&& n2.geom and st_3dintersects(n1.geom, n2.geom) and n1.id<1000;
COPY nuclei1(geom) FROM '/home/teng/project/HiSPEED/src/wkt_n_nv10_nu5000_s10.dt' WITH DELIMITER AS '|';

select pg_relation_size('actor');

select pg_size_pretty(pg_relation_size('nuclei2'));


CREATE TABLE nuclei1(id SERIAL PRIMARY KEY, geom geometry);
CREATE TABLE nuclei2(id SERIAL PRIMARY KEY, geom geometry);
CREATE TABLE vessel(id SERIAL PRIMARY KEY, geom geometry);
create index global_nuclei1_gix on nuclei1 using GIST(geom);
create index global_nuclei2_gix on nuclei2 using GIST(geom);
create index global_vessel_gix on vessel using GIST(geom);
COPY vessel(geom) FROM '/home/teng/project/HiSPEED/src/wkt_v_nv200_nu200_s10.wkt' WITH DELIMITER AS '|';
COPY nuclei1(geom) FROM '/home/teng/project/HiSPEED/src/wkt_n_nv200_nu200_s10.wkt' WITH DELIMITER AS '|';
COPY nuclei2(geom) FROM '/home/teng/project/HiSPEED/src/wkt_n2_nv200_nu200_s10.wkt' WITH DELIMITER AS '|';


select	g1.id,g2.id,ST_3DDistance(g1.geom, g2.geom) as dist 
from	nuclei1 g1, nuclei2 g2 
where	st_expand(g1.geom, 10) &&& n2.geom and n1.id=1
order by dist asc limit 1;