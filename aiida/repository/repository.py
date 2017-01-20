# -*- coding: utf-8 -*-

from abc import ABCMeta, abstractmethod

class Repository(object):
    """
    """
    __metaclass__ = ABCMeta

    @abstractmethod
    def __init__(self, repo_config):
        """
        :param repo_config: dictionary with configuration details for repository
        """
        raise NotImplementedError

    @abstractmethod
    def get_name(self):
        """
        Return the name of the repository which is a human-readable label

        :return name: the human readable label associated with this repository
        """
        raise NotImplementedError

    @abstractmethod
    def get_uuid(self):
        """
        Return the UUID identifying the repository. How and where it is stored
        is implementation specific
        :return: string with repository UUID
        """
        raise NotImplementedError

    @abstractmethod
    def put_object(self, key, content, stop_if_existing=False):
        """
        Add the object to the repository with the specified content if 
        it does not exist. Overwrite if exists and stop_if_existing is False.
        Raise an exception if stop_if_existing is True and the object already exists.

        :param key: fully qualified identifier for the object within the repository
        :param content: file or filelike
        :param stop_if_existing:
        """
        raise NotImplementedError

    @abstractmethod
    def put_new_object(self, source, node_uuid, relpath):
        """
        Add the object to the repository with some autogenerated name.
        Can (optionally) use the information from node_uuid (string) and relpath (string),
        but it is not needed.
        Must return the key used to store the object.
        It never overwrites existing data. If you want to overwrite an existing object,
        call put_object passing the correct key.

        :return: the key of the newly generated object
        """
        raise NotImplementedError

    @abstractmethod
    def exists(self, key):
        """
        Return True or False. In swift, use head_object() and catch 404 exception.
        """
        raise NotImplementedError

    @abstractmethod
    def get_object(self, key):
        """
        Get the content
        """
        raise NotImplementedError

    @abstractmethod
    def del_object(self, key):
        """
        Delete the object. This class is still valid, if you put again you recreate the object.
        """
        raise NotImplementedError