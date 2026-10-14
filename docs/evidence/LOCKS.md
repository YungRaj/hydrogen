# Evidence locks

Locks prove that a preregistered analysis was not tuned after its outcomes
were seen. They do **not** freeze the working tree.

## Completed campaigns
`docs/evidence/*_lock.json` files record the source hashes an analysis used.
[`lock_provenance.json`](lock_provenance.json) anchors each one to the git
commit where those hashes verify, and
`test_historical_locks_verify_at_their_recorded_commits` re-checks them
against that commit with `git show`. Code at HEAD may change freely.

To re-run a frozen analysis, check out its recorded commit. The NEB runner
refuses a modified implementation unless you pass `--exploratory`; such runs
are labelled `exploratory_modified_implementation` and can never pass
acceptance.

Never edit a completed lock file: result artifacts record its SHA-256.

## Future campaigns: lock at launch
`prospective_policy_lock_current.json` holds the campaign design but no source
hash until launch. Immediately before the first `prereg-` batch run:

```bash
python run_divide_conquer_pilot.py lock-policy
```

This records the digest of the runner's **full import closure**
(`pipeline.evidence.provenance.closure_digest`), the commit and the time.
Every `prepare` and the confirmatory analysis then require an identical
closure. Re-locking is refused once confirmatory batches exist (`--force`
overrides, deliberately splitting the campaign). After the campaign
completes, add its lock and commit to `lock_provenance.json`.
