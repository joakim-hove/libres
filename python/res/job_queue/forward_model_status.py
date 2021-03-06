#  Copyright (C) 2018  Statoil ASA, Norway.
#
#  The file 'forward_model_status.py' is part of ERT - Ensemble based Reservoir Tool.
#
#  ERT is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  ERT is distributed in the hope that it will be useful, but WITHOUT ANY
#  WARRANTY; without even the implied warranty of MERCHANTABILITY or
#  FITNESS FOR A PARTICULAR PURPOSE.
#
#  See the GNU General Public License at <http://www.gnu.org/licenses/gpl.html>
#  for more details.
import os.path
import json
import datetime
import time
import sys

def _serialize_date(dt):
    if dt is None:
        return None

    return time.mktime(dt.timetuple())


def _deserialize_date(serial_dt):
    if serial_dt is None:
        return None

    time_struct = time.localtime(serial_dt)
    return datetime.datetime(*time_struct[0:6])


class ForwardModelJobStatus(object):

    def __init__(self, name, start_time = None, end_time = None, status = "Waiting", error=None):
        self.start_time = start_time
        self.end_time = end_time
        self.name = name
        self.status = status
        self.error = error


    @classmethod
    def load(cls,data):
        start_time = _deserialize_date(data["start_time"])
        end_time = _deserialize_date(data["end_time"])
        name = data["name"]
        status = data["status"]
        error = data["error"]

        return cls(name,
                   start_time=start_time,
                   end_time=end_time,
                   status=status,
                   error=error)


    def __str__(self):
        return "name:{} start_time:{}  end_time:{}  status:{}  error:{} ".format(self.name, self.start_time, self.end_time, self.status, self.error)


    def dump_data(self):
        return {"name" : self.name,
                "status" : self.status,
                "error" : self.error,
                "start_time" : _serialize_date(self.start_time),
                "end_time" : _serialize_date(self.end_time)}

class ForwardModelStatus(object):
    STATUS_FILE = "status.json"


    def __init__(self, run_id, start_time, end_time = None):
        self.run_id = run_id
        self.start_time = start_time
        self.end_time = end_time
        self._jobs = []

    @classmethod
    def try_load(cls, status_file):
        fp = open(status_file)
        data = json.load(fp)

        start_time = _deserialize_date(data["start_time"])
        end_time = _deserialize_date(data["end_time"])
        status = cls(data["run_id"],
                     start_time,
                     end_time=end_time)

        for job in data["jobs"]:
            status.add_job(ForwardModelJobStatus.load(job))

        return status


    @classmethod
    def load(cls, path, num_retry=10):
        sleep_time = 0.10
        attempt = 0
        status_file = os.path.join(path, cls.STATUS_FILE)
        while attempt < num_retry:
            try:
                status = cls.try_load(status_file)
                return status
            except:
                attempt += 1
                if attempt < num_retry:
                    time.sleep(sleep_time)

        return None

    @property
    def jobs(self):
        return self._jobs


    def add_job(self, job):
        self._jobs.append(job)


    def dump(self, filename = None):
        if filename is None:
            status_file = self.STATUS_FILE
        else:
            status_file = filename

        data = {"run_id" : self.run_id,
                "start_time" : _serialize_date(self.start_time),
                "end_time" : _serialize_date(self.end_time)}
        jobs = []
        for job in self.jobs:
            jobs.append( job.dump_data() )

        data["jobs"] = jobs
        with open(status_file, "w") as fp:
            json.dump(data, fp)


    def complete(self):
        self.end_time = datetime.datetime.now()
        self.dump( )
