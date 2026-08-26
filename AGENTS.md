# Codex Project Instructions

## General Objective

Work on this repository carefully, efficiently, and with minimal unnecessary token usage.

Preserve existing functionality unless the current task explicitly requires changing it.

Before making changes, inspect the relevant existing implementation and understand how the affected components interact. Do not assume architecture, file locations, APIs, models, database structures, or behavior without checking the repository.

## Cost-Efficient Agent Strategy

Use the currently selected parent model as the primary implementation and integration agent.

When multi-agent delegation and model selection are available:

- Prefer GPT-5.6 Luna for narrow, well-defined, low-risk, mechanical work.
- Prefer GPT-5.6 Terra for architecture, implementation, debugging, integration, and decisions requiring substantial reasoning.
- Do not use a higher-cost model when a lower-cost model can reliably complete the work.

Good tasks for Luna include:

- repository searches
- finding references or usages
- localization/string audits
- repetitive edits
- straightforward renames
- formatting
- simple isolated tests
- checking for missing cases
- identifying files that contain a particular pattern
- mechanical migrations
- basic documentation updates
- simple validation work

Keep Terra responsible for:

- architecture
- cross-file design decisions
- state management
- data flow
- backend/frontend coordination
- database or schema changes
- ambiguous bugs
- concurrency
- authentication or permissions
- complicated SwiftUI behavior
- complicated API behavior
- integration
- reviewing delegated changes
- final verification

When delegating:

- Give each subagent only the context necessary for its specific task.
- Avoid sending the entire project history or conversation to a subagent unnecessarily.
- Keep subagent assignments narrow and independently verifiable.
- Parallelize independent searches or audits when useful.
- Review delegated work before incorporating it.
- The parent agent remains responsible for the final result.

Do not delegate merely for the sake of delegation. Small tasks that are faster to perform directly should be performed directly.

## Workflow

For non-trivial tasks:

1. Inspect the relevant repository structure and existing implementation.
2. Determine the actual root cause or required architecture before editing.
3. Identify all affected surfaces.
4. Make the smallest coherent set of changes necessary.
5. Check for related implementations that must remain consistent.
6. Run relevant tests, builds, type checks, linters, or validation.
7. Fix problems caused by the changes.
8. Review the final diff for unintended modifications.
9. Report what changed and how it was verified.

Do not stop after merely writing code when the repository provides a reasonable way to validate it.

## Existing Functionality

Do not:

- remove existing functionality unless explicitly requested
- perform unrelated refactors
- rename unrelated files or APIs
- change dependencies without a reason
- replace working implementations simply because another approach is preferred
- silently change user-facing behavior outside the requested scope

When modifying an existing feature, inspect related implementations so behavior remains consistent across the application.

## Repository Searches

When a task could affect multiple implementations, search the repository rather than assuming the first occurrence is the only occurrence.

Examples include:

- localization strings
- navigation destinations
- shared models
- API endpoints
- duplicated mobile/web behavior
- feature flags
- database fields
- analytics
- notifications
- tests
- accessibility labels

## Testing and Validation

Use the repository's existing testing and build systems.

Prefer targeted validation during development, followed by the broadest reasonable validation before completion.

If a test or build cannot be run, clearly state:

- what could not be run
- why
- what was validated instead

Never claim that something passed unless it was actually run successfully.

## Errors Discovered During Work

If you encounter an existing error directly related to the area being modified, investigate and fix it when doing so is safe and within scope.

Do not ignore compilation errors, broken references, failing tests, or obvious inconsistencies caused by the current changes.

Do not expand the task into unrelated cleanup.

## Final Review

Before declaring the task complete:

- review all changed files
- make sure there are no accidental changes
- confirm requested behavior is implemented
- check relevant existing functionality
- run available validation
- resolve issues introduced by the work

Provide a concise final summary containing:

1. What changed
2. Important implementation decisions
3. Validation/tests performed
4. Anything that still requires attention
