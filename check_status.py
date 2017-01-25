#!/usr/bin/env python3

import click
from sqlalchemy import and_, asc, desc, func, MetaData, Table, or_
from sqlalchemy.sql import select
import yaml

from htsohm.db.__init__ import engine

meta = MetaData(bind=engine)
materials = Table('materials', meta, autoload=True)
mutation_strengths = Table('mutation_strengths', meta, autoload=True)

@click.group()
def data():
    pass

@data.command()
def run_ids():
    cols = [materials.c.run_id,
            func.max(materials.c.generation),
            func.count(materials.c.id)]
    rows = or_(materials.c.retest_passed == None, materials.c.retest_passed == True)
    sort = materials.c.run_id
    s = select(cols, rows).group_by(sort).order_by(asc(sort))
    print('\nrun-id\t\t\t\tgenerations\tmaterial')
    result = engine.execute(s)
    for row in result:
        print('%s\t%s\t\t%s' % (row[0], row[1], row[2]))
    result.close()

@data.command()
@click.argument('run_id')
@click.argument('parameter')
def maximum(run_id, parameter):
    param = getattr(materials.c, parameter)
    cols = [materials.c.uuid, param]
    rows = and_(materials.c.run_id == run_id,
            or_(materials.c.retest_passed == None,
                materials.c.retest_passed == True))
    s = select(cols, rows).order_by(desc(param)).limit(1)
    print('\nuuid\t\t\t\t\t%s' % parameter)
    result = engine.execute(s)
    for row in result:
        print('%s\t%s' % (row[0], row[1]))
    result.close()

@data.command()
@click.argument('uuid')
def material(uuid):
    cols = [materials.c.uuid,
            materials.c.parent_id,
            materials.c.ga_absolute_volumetric_loading,
            materials.c.sa_volumetric_surface_area,
            materials.c.vf_helium_void_fraction,
            materials.c.generation,
            materials.c.run_id]
    rows = and_(materials.c.uuid == uuid,
            or_(materials.c.retest_passed == None,
                materials.c.retest_passed == True))
    s = select(cols, rows)
    print(
        '\nuuid\t\t\t\t\tparent\tgas adsorption (cc/cc)\tsurface area (m2/cc)' +
        '\tvoid fraction\tgeneration\trun'
    )
    result = engine.execute(s)
    for row in result:
        print(
            '%s\t%s\t%s\t\t%s\t\t\t' %  (row[0], row[1], row[2], row[3]) +
            '%s\t%s\t\t%s' %  (row[4], row[5], row[6])    
        )
    result.close()

@data.command()
@click.argument('uuid')
def find_children(uuid):
    cols = [materials.c.id]
    rows = [materials.c.uuid == uuid]
    result = engine.execute(select(cols, *rows))
    for row in result:
        parent_id = row[0]
    result.close()

    cols = [materials.c.uuid]
    rows = and_(materials.c.parent_id == parent_id,
            or_(materials.c.retest_passed == None,
                materials.c.retest_passed == True))
    print('\nchildren of %s :' % uuid)
    result = engine.execute(select(cols, rows))
    for row in result:
        print('\t%s' % row[0])
    result.close()

if __name__ == '__main__':
    data()
