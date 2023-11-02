from string import ascii_lowercase as LETTERS
from typing import Optional

BIN_SIZE = 100_000


class Labeller:
    def __init__(self):
        self.chr: Optional[str] = None
        self.bin: Optional[int] = None
        self.cnt: int = 0

    def _cnt_to_ascii(self) -> str:
        result, cnt = "", self.cnt
        while cnt > 0:
            cnt, rem = divmod(cnt, len(LETTERS))
            result += LETTERS[rem]
        return result

    def next_label(self, chrom: str, start: int) -> str:
        curbin = start // BIN_SIZE
        if (self.chr == chrom) and (self.bin == curbin):
            self.cnt += 1
            label = self._cnt_to_ascii()
        else:
            label = 'a'
            self.chr = chrom
            self.bin = curbin
            self.cnt = 0
        return f'X_{chrom.strip("chr")}_{curbin}_{label}'
