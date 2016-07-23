from __future__ import division

import os, time
import libtbx.load_env
from libtbx.utils import Sorry

from xfel.ui.db.trial import Trial
from xfel.ui.db.run import Run
from xfel.ui.db.rungroup import Rungroup
from xfel.ui.db.tag import Tag
from xfel.ui.db.job import Job
from xfel.ui.db.stats import Stats
from xfel.ui.db.experiment import Cell, Bin, Isoform, Event

from xfel.ui.db import get_db_connection

try:
  from MySQLdb import OperationalError
except ImportError:
  raise Sorry('Mysql not available')

from xfel.command_line.experiment_manager import initialize as initialize_base

class initialize(initialize_base):
  expected_tables = ["run", "job", "rungroup", "trial", "tag", "run_tag", "event", "trial_rungroup",
                     "imageset", "imageset_event", "beam", "detector", "experiment",
                     "crystal", "cell", "cell_bin", "bin", "isoform"]

  def __init__(self, params, dbobj):
    initialize_base.__init__(self, params, dbobj, interactive = False, drop_tables = None)

  def create_tables(self, sql_path = None):
    if sql_path is None:
      sql_path = os.path.join(libtbx.env.find_in_repositories("xfel/ui/db"), "schema.sql")

    return initialize_base.create_tables(self, sql_path)

  def set_up_columns_dict(self, app):
    columns_dict = {}
    for table in self.expected_tables:
      table_name = "%s_%s" % (self.params.experiment_tag, table)
      query = "SHOW COLUMNS FROM `%s`" % (table_name)
      cursor = app.execute_query(query)
      columns_dict[table_name] = [c[0] for c in cursor.fetchall() if c[0] != 'id']
    return columns_dict

class xfel_db_application(object):
  def __init__(self, params, drop_tables = False, verify_tables = False):
    self.params = params
    self.query_count = 0
    dbobj = get_db_connection(params)
    self.init_tables = initialize(params, dbobj) # only place where a connection is held

    if drop_tables:
      self.drop_tables()

    if verify_tables and not self.verify_tables():
      self.create_tables()
      print 'Creating experiment tables...'
      if not self.verify_tables():
        raise Sorry("Couldn't create experiment tables")

    self.columns_dict = self.init_tables.set_up_columns_dict(self)

  def execute_query(self, query, commit = False, verbose = False):
    if verbose:
      from time import time
      st = time()
      self.query_count += 1

    retry_count = 0
    retry_max = 10
    sleep_time = 0.1
    while retry_count < retry_max:
      try:
        dbobj = get_db_connection(self.params)
        cursor = dbobj.cursor()
        cursor.execute(query)
        if commit:
          dbobj.commit()

        if verbose:
          print 'Query % 6d SQLTime Taken = % 10.6f seconds' % (self.query_count, time() - st), query[:min(len(query),145)]
        return cursor
      except OperationalError, e:
        if "Can't connect to MySQL server" not in str(e):
          raise e
        retry_count += 1
        print "Couldn't connect to MYSQL, retry", retry_count
        time.sleep(sleep_time)
        sleep_time *= 2
      except Exception, e:
        print "Couldn't execute MYSQL query.  Query:"
        print query
        print "Exception:"
        print str(e)
        raise e
    raise Sorry("Couldn't execute MYSQL query. Too many reconnects. Query: %s"%query)

  def list_lcls_runs(self):
    from xfel.xpp.simulate import file_table
    query = "https://pswww.slac.stanford.edu/ws-auth/dataexport/placed?exp_name=%s" % (self.params.experiment)
    FT = file_table(self.params, query, enforce80=self.params.web.enforce80, enforce81=self.params.web.enforce81)
    return FT.get_runs()

  def verify_tables(self):
    return self.init_tables.verify_tables()

  def create_tables(self):
    return self.init_tables.create_tables()

  def drop_tables(self):
    return self.init_tables.drop_tables()

  def create_trial(self, d_min = 1.5, n_bins = 10, **kwargs):
    # d_min and n_bins only used if isoforms are in this trial

    trial = Trial(self, **kwargs)
    if trial.target_phil_str is not None:
      from iotbx.phil import parse
      backend = ['labelit', 'dials'][['cxi.xtc_process', 'cctbx.xfel.xtc_process'].index(self.params.dispatcher)]
      if backend == 'labelit':
        from spotfinder.applications.xfel import cxi_phil
        trial_params = cxi_phil.cxi_versioned_extract().persist.phil_scope.fetch(parse(trial.target_phil_str)).extract()
        isoforms = trial_params.isoforms
      elif backend == 'dials':
        from xfel.command_line.xtc_process import phil_scope
        trial_params = phil_scope.fetch(parse(trial.target_phil_str)).extract()
        isoforms = trial_params.indexing.stills.isoforms
      else:
        assert False
      if len(isoforms) > 0:
        for isoform in isoforms:
          print "Creating isoform", isoform.name
          db_isoform = Isoform(self,
                               name = isoform.name,
                               trial_id = trial.id)
          a, b, c, alpha, beta, gamma = isoform.cell.parameters()
          cell = self.create_cell(cell_a = a,
                                  cell_b = b,
                                  cell_c = c,
                                  cell_alpha = alpha,
                                  cell_beta = beta,
                                  cell_gamma = gamma,
                                  lookup_symbol = isoform.lookup_symbol,
                                  isoform_id = db_isoform.id)
          from cctbx.crystal import symmetry

          cs = symmetry(unit_cell = isoform.cell,space_group_symbol=str(isoform.lookup_symbol))
          mset = cs.build_miller_set(anomalous_flag=False, d_min=d_min)
          binner = mset.setup_binner(n_bins=n_bins)
          for i in binner.range_used():
            d_max, d_min = binner.bin_d_range(i)
            Bin(self,
                number = i,
                d_min = d_min,
                d_max = d_max,
                total_hkl = binner.counts_complete()[i],
                cell_id = cell.id)
    return trial

  def get_trial_isoforms(self, trial_id):
    where = "WHERE trial_id = %d"%trial_id
    return self.get_all_x(Isoform, "isoform", where)

  def create_cell(self, **kwargs):
    return Cell(self, **kwargs)

  def get_cell(self, cell_id = None, name = None):
    assert [cell_id, name].count(None) == 1
    if name is not None:
      query = "SELECT id FROM `%s_cell` WHERE name = '%s'"%(self.params.experiment_tag, name)
      cursor = self.execute_query(query)
      results = cursor.fetchall()
      assert len(results) in [0,1]
      if len(results) == 0:
        return None
      cell_id = int(results[0][0])

    return Cell(self, cell_id=cell_id)

  def get_cell_bins(self, cell_id):
    query = "SELECT id FROM `%s_bin` WHERE cell_id = %d" % \
            (self.params.experiment_tag, cell_id)
    cursor = self.execute_query(query)
    ids = [str(i[0]) for i in cursor.fetchall()]
    if len(ids) == 0:
      return []
    where = "WHERE id IN (%s)" % ", ".join(ids)
    return self.get_all_x(Bin, 'bin', where)

  def get_all_x(self, cls, name, where = None):
    table_name = "%s_%s" % (self.params.experiment_tag, name)
    columns = ["%s.%s"%(name, c) for c in self.columns_dict[table_name]]
    query = "SELECT %s.id, %s FROM `%s` %s" % (name, ", ".join(columns), table_name, name)
    if where is not None:
      query += " " + where
    cursor = self.execute_query(query)
    results = []
    for row in cursor.fetchall():
      d = {key:value for key, value in zip(self.columns_dict[table_name], row[1:])}
      d["%s_id"%name] = row[0]
      results.append(cls(self, **d))
    return results

  def get_all_x_with_subitems(self, cls, name, where = None, sub_items = None):
    """ Assemble a list of db_proxy objects, where each one references a sub object
        @param cls principal class
        @param name table name not including experiment tag
        @param where optional constraints on what to get
        @param sub_items: array of tuples (cls,name) of items this table references
        @return list of assembled db_proxy objects"""

    # Assemble a list of all columns to be extracted via a big join. Then entries in the
    # columns array will be of the form "name.column", IE: tag.comment
    table_name = "%s_%s" % (self.params.experiment_tag, name)
    columns = ["%s.%s"%(name, c) for c in self.columns_dict[table_name]]
    columns.append("%s.id"%name)

    if sub_items is None:
      sub_items = []
    sub_table_names = ["%s_%s"%(self.params.experiment_tag, i[1]) for i in sub_items]
    for i, sub_item in enumerate(sub_items):
      scls, sname = sub_item
      columns.extend(["%s.%s"%(sname, c) for c in self.columns_dict[sub_table_names[i]]])
      columns.append("%s.id"%sname)

    # the main item being extracted is in the FROM statement and is given a nickname which
    # is the table name without the experiment tag
    query = "SELECT %s FROM `%s` %s" % (", ".join(columns), table_name, name)

    # Join statements to bring in the sub tables
    for i, sub_item in enumerate(sub_items):
      scls, sname = sub_item
      query += " JOIN `%s` %s ON %s.id = %s.%s_id"% (
        sub_table_names[i], sname, sname, name, sname)

    if where is not None:
      query += " " + where
    cursor = self.execute_query(query)

    results = []
    for row in cursor.fetchall():
      # Each row will be a complete item and sub items in column form. Assemble one
      # dictionary (d) for the main item and a dictionary of dictionaries (sub_ds)
      # for each of the sub items
      d = {}
      sub_ds = {sub_item[1]:(sub_item[0], {}) for sub_item in sub_items}
      for key, value in zip(columns, row):
        n, c = key.split('.') # nickname n, column name c
        if n == name:
          d[c] = value # this column came from the main table
        else:
          sub_ds[n][1][c] = value # this column came from a sub table

      # pop the id column as it is passed as name_id to the db_proxy class (ie Job(job_id = 2))
      _id = d.pop("id")
      d["%s_id"%name] = _id
      results.append(cls(self, **d)) # instantiate the main class
      for sub_d_n, sub_d in sub_ds.iteritems():
        _id = sub_d[1].pop("id")
        sub_d[1]["%s_id"%sub_d_n] = _id
        setattr(results[-1], sub_d_n, sub_d[0](self, **sub_d[1])) # instantiate the sub items
    return results

  def get_trial(self, trial_id = None, trial_number = None):
    assert [trial_id, trial_number].count(None) == 1
    if trial_id is None:
      trials = [t for t in self.get_all_trials() if t.trial == trial_number]
      assert len(trials) == 1
      return trials[0]
    else:
      return Trial(self, trial_id)

  def get_trial_rungroups(self, trial_id, only_active = False):
    query = "SELECT rungroup_id FROM `%s_trial_rungroup` WHERE `%s_trial_rungroup`.trial_id = %d" % \
            (self.params.experiment_tag, self.params.experiment_tag, trial_id)
    cursor = self.execute_query(query)
    rungroups = [Rungroup(self, i[0]) for i in cursor.fetchall()]
    if only_active:
      return [rg for rg in rungroups if rg.active]
    else:
      return rungroups

  def get_trial_runs(self, trial_id):
    rungroups = self.get_trial_rungroups(trial_id, only_active=True)
    runs = []
    run_ids = []
    for rungroup in rungroups:
      for run in rungroup.runs:
        if run.id not in run_ids:
          run_ids.append(run.id)
          runs.append(run)
    return runs

  def get_all_trials(self, only_active = False):
    if only_active:
      return [t for t in self.get_all_x(Trial, "trial") if t.active]
    else:
      return self.get_all_x(Trial, "trial")

  def create_run(self, **kwargs):
    return Run(self, **kwargs)

  def get_run(self, run_id = None, run_number = None):
    assert [run_id, run_number].count(None) == 1
    if run_id is not None:
      return Run(self, run_id)
    runs = [r for r in self.get_all_runs() if r.run == run_number]
    if len(runs) == 0:
      raise Sorry("Couldn't find run %d"%run_number)
    assert len(runs) == 1
    return runs[0]

  def get_all_runs(self):
    return self.get_all_x(Run, "run")

  def create_rungroup(self, **kwargs):
    return Rungroup(self, **kwargs)

  def get_rungroup(self, rungroup_id):
    return Rungroup(self, rungroup_id)

  def get_all_rungroups(self, only_active = False):
    if only_active:
      return [rg for rg in self.get_all_x(Rungroup, "rungroup") if rg.active]
    else:
      return self.get_all_x(Rungroup, "rungroup")

  def create_tag(self, **kwargs):
    return Tag(self, **kwargs)

  def get_tag(self, tag_id):
    return Tag(self, tag_id)

  def get_run_tags(self, run_id):
    query = "SELECT tag_id from `%s_run_tag` WHERE `%s_run_tag`.run_id = %d" % \
            (self.params.experiment_tag, self.params.experiment_tag, run_id)
    cursor = self.execute_query(query)
    tag_ids = [str(i[0]) for i in cursor.fetchall()]
    if len(tag_ids) == 0:
      return []
    where = "WHERE id IN (%s)" % ", ".join(tag_ids)
    return self.get_all_x(Tag, 'tag', where)

  def get_all_tags(self):
    return self.get_all_x(Tag, "tag")

  def delete_x(self, item, item_id):
    query = "DELETE FROM `%s` WHERE id = %d"%(item.table_name, item_id)
    self.execute_query(query, commit = True)

  def delete_tag(self, tag = None, tag_id = None):
    assert [tag, tag_id].count(None) == 1

    if tag_id is None:
      tag_id = tag.tag_id
    else:
      tag = self.get_tag(tag_id)

    query = "DELETE FROM `%s_run_tag` WHERE tag_id = %d" % (self.params.experiment_tag, tag_id)
    self.execute_query(query, commit=True)

    self.delete_x(tag, tag_id)

  def create_job(self, **kwargs):
    return Job(self, **kwargs)

  def get_job(self, job_id):
    return Job(self, job_id)

  def get_all_jobs(self):
    return self.get_all_x_with_subitems(Job, "job", sub_items = [(Trial, 'trial'), (Run, 'run'), (Rungroup, 'rungroup')])

  def delete_job(self, job = None, job_id = None):
    assert [job, job_id].count(None) == 1
    if job_id is None:
      job_id = job.job_id
    else:
      job = self.get_job(job_id)

    self.delete_x(job, job_id)

  def get_all_events(self, trial = None, runs = None, only_indexed = True, isoform = None, where = None):
    tag = self.params.experiment_tag
    if where is None:
      final_where = ""
    else:
      final_where = where.strip()
    where = ""
    if only_indexed or isoform is not None:
      where = " JOIN `%s_imageset_event` is_e ON event.id = is_e.event_id"%tag
    if trial is not None:
      if runs is None:
        runs = trial.runs
      if len(runs) == 0:
        return []
      if isoform is None:
        where += " WHERE"
      else:
        where += """ JOIN `%s_imageset` imgset ON imgset.id = is_e.imageset_id
                     JOIN `%s_experiment` exp ON exp.imageset_id = imgset.id
                     JOIN `%s_crystal` crystal ON crystal.id = exp.crystal_id
                     JOIN `%s_cell` cell ON cell.id = crystal.cell_id
                     JOIN `%s_isoform` isoform ON isoform.id = cell.isoform_id
                     WHERE isoform.id = %s AND"""%(tag, tag, tag, tag, tag, isoform.id)
      where += " event.trial_id = %d AND event.run_id in (%s)" % (
        trial.id, ", ".join([str(r.id) for r in runs]))

      if 'rungroup_id' in self.columns_dict["%s_%s" % (tag, 'event')]: # some backwards compatibility, as event.rungroup_id was added late to the schema
        rungroups = ", ".join([str(rg.id) for rg in trial.rungroups])
        if len(rungroups) > 0:
          where += " AND event.rungroup_id in (%s)"%rungroups

    if len(where) > 0 and len(final_where) > 0:
      final_where = "AND " + final_where.lstrip("WHERE")
    return self.get_all_x(Event, "event", where + " " + final_where)

  def get_stats(self, **kwargs):
    return Stats(self, **kwargs)
