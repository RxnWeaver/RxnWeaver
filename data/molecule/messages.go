package molecule

// RequestType enumerates the requests understood by a molecule.
type RequestType uint8

// StatusType enumerates the response statuses sent by a molecule.
type StatusType uint16

// InMessage is a message sent to a molecule by an external agent.
//
// An in-message comprises a request, an optional cookie value and a
// dynamic payload.  The meaning of the cookie depends on the request,
// but usually represents a request ID.  Similarly, the actual
// contents of the payload depend on the request.
//
// Each request also includes the channel over which the result
// message has to be transmitted.  This enables the molecules to
// operate in a stateless manner, as far as the communication is
// concerned.
//
// Thus, it is highly imperative that other agents that correspond
// with a molecule be aware of what requests molecules understand, and
// what payloads are to be delivered as part of the message.
type InMessage struct {
	Request    RequestType
	Cookie     uint64
	OutChannel chan OutMessage
	Payload    interface{}
}

// OutMessage is a message sent by a molecule in response to an
// in-message.
//
// An out-message comprises a status (result code), an optional cookie
// value and a dynamic payload.  The meaning of the cookie depends on
// the request, but usually is the same as the corresponding request
// ID.  Similarly, the actual contents of the payload depend on the
// request.
//
// Thus, it is highly imperative that other agents that correspond
// with a molecule be aware of what responses molecules send, and what
// payloads are delivered as part of the message.
type OutMessage struct {
	Status  StatusType
	Cookie  uint64
	Payload interface{}
}

// Constants representing the requests understood by a molecule.
const (
	ReqNone RequestType = iota // Do not use this.
	ReqAddAtom
	ReqAddBond
	ReqSetAtomAttribute
	ReqAddTag
)

// Constants representing the outcome status of a request processed by
// a molecule.
const (
	StSuccess StatusType = iota
	StNotFound
	StAlreadyExists
	StIncorrectParameter
)
