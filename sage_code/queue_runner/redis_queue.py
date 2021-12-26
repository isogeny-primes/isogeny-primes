import json
from pathlib import Path
from typing import Optional, Union

import redis
from pydantic import BaseSettings, Field
from .sage_json_converter import sage_converter


class RedisSettings(BaseSettings):
    host: str = Field("localhost", env="REDIS_HOST")
    port: int = Field(6379, env="REDIS_PORT")
    db: int = Field(0, env="REDIS_DB")
    password: Optional[str] = Field(None, env="REDIS_PASSWORD")
    socket_timeout: Optional[int] = Field(None, env="REDIS_SOCKET_TIMEOUT")
    socket_connect_timeout: Optional[int] = Field(
        None, env="REDIS_SOCKET_CONNECT_TIMEOUT"
    )
    socket_keepalive: Optional[int] = Field(None, env="REDIS_SOCKET_KEEPALIVE")
    socket_keepalive_options: Optional[str] = Field(
        None, env="REDIS_SOCKET_KEEPALIVE_OPTIONS"
    )
    connection_pool: Optional[str] = Field(None, env="REDIS_CONNECTION_POOL")
    unix_socket_path: Optional[str] = Field(None, env="REDIS_UNIX_SOCKET_PATH")
    encoding: str = Field("utf-8", env="REDIS_ENCODING")
    encoding_errors: str = Field("strict", env="REDIS_ENCODING_ERRORS")
    charset: Optional[str] = Field(None, env="REDIS_CHARSET")
    errors: Optional[str] = Field(None, env="REDIS_ERRORS")
    decode_responses: bool = Field(False, env="REDIS_DECODE_RESPONSES")
    retry_on_timeout: bool = Field(False, env="REDIS_RETRY_ON_TIMEOUT")
    ssl: bool = Field(False, env="REDIS_SSL")
    ssl_keyfile: Optional[str] = Field(None, env="REDIS_SSL_KEYFILE")
    ssl_certfile: Optional[str] = Field(None, env="REDIS_SSL_CERTFILE")
    ssl_cert_reqs: str = Field("required", env="REDIS_SSL_CERT_REQS")
    ssl_ca_certs: Optional[str] = Field(None, env="REDIS_SSL_CA_CERTS")
    ssl_check_hostname: bool = Field(False, env="REDIS_SSL_CHECK_HOSTNAME")
    max_connections: Optional[int] = Field(None, env="REDIS_MAX_CONNECTIONS")
    single_connection_client: bool = Field(False, env="REDIS_SINGLE_CONNECTION_CLIENT")
    health_check_interval: int = Field(0, env="REDIS_HEALTH_CHECK_INTERVAL")
    client_name: Optional[str] = Field(None, env="REDIS_CLIENT_NAME")
    username: Optional[str] = Field(None, env="REDIS_USERNAME")


class BaseFifoQueue:
    def _push(self, name, *values, block_size: int = 1000):
        raise NotImplementedError

    def _pop(self, name, count=None):
        raise NotImplementedError

    def _range(self, name, start, stop):
        raise NotImplementedError

    def push(self, name, *values):
        values = [json.dumps(v, default=sage_converter) for v in values]
        self._push(name, *values)

    def push_file(self, name, file_name):
        with open(file_name, "r") as file:
            data = file.read()
        self.push(name, *json.loads(data))

    def pop(self, name):
        value = self._pop(name)
        if value is None:
            raise IndexError("pop from empty queue")
        return json.loads(value.decode())

    def clear(self, name):
        raise NotImplementedError

    def len(self, name):
        raise NotImplementedError

    def range(self, name, start, stop):
        result = self._range(name, start, stop)
        return [json.loads(r.decode()) for r in result]

    def dump(self, name, file_name, block_size: int = 1000):
        result = []
        for block in range(0, self.len(name), block_size):
            result.extend(self.range(name, block, block + block_size))
        result_str = ",\n".join(
            json.dumps(r, separators=(",", ":"), sort_keys=True) for r in result
        )
        result_str = f"[\n{result_str}\n]"
        with open(file_name, "w") as file:
            file.write(result_str)


class RedisFifoQueue(BaseFifoQueue):
    def __init__(
        self, redis_settings: Optional[Union[RedisSettings, str, Path]] = None
    ):
        if not redis_settings:
            redis_settings = RedisSettings()

        if isinstance(redis_settings, (str, Path)):
            redis_settings = RedisSettings(_env_file=str)

        self._redis: redis.Redis = redis.Redis(**redis_settings.dict())

    def _push(self, name, *values, block_size: int = 1000):
        if not block_size:
            self._redis.rpush(name, *values)
            return

        for block in range(0, len(values), block_size):
            self._redis.rpush(name, *values[block : block + block_size])

    def _pop(self, name, count=None):
        return self._redis.lpop(name, count)

    def _range(self, name, start, stop):
        return self._redis.lrange(name, start, stop - 1)

    def clear(self, name):
        self._redis.delete(name)

    def len(self, name):
        return self._redis.llen(name)


if __name__ == "__main__":
    queue = RedisFifoQueue()
    queue.push("test", {"a": "value"})
    print(queue.pop("test"))
    print(queue.pop("test"))
